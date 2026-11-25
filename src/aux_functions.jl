#--------------------------------------------------------------#
#-Networks, hurricane regimes and the interaction matrix-------#
#--------------------------------------------------------------#

# Functions shared by the scripts in scripts/03_simulations: loading one
# multilayer network, the damage fraction of each hurricane category, the
# empirical hurricane regime of each island and the interaction matrix of
# the population dynamics.
#
# Inputs:  the two incidence matrices of a network, read as DataFrames
#          src/empirical_lambdas.csv
#          src/empirical_pmfs.csv
# Outputs: none. The functions return the inputs of the simulations.

using CSV
using DataFrames
using LinearAlgebra

#--------------------------------------------------------------#
#-Species and binary layers of one network---------------------#
#--------------------------------------------------------------#

# Read one multilayer network from its pollination (B_df) and seed-dispersal
# (S_df) incidence matrices, with plants in the first column and one column
# per bird. Every positive entry is an interaction. Keep only the species
# with at least one interaction in either layer, plants first and then
# birds. A species in both layers is "Multilayer" (for birds, the dual-role
# birds). The others are "Pollination" or "Seed-dispersal".
# Returns a NamedTuple with network, species, category and the binary
# matrices B and S of the species kept.
function load_network(; B_df::DataFrame, S_df::DataFrame, network::String)
    plants = String.(B_df[!, 1])
    birds = String.(names(B_df)[2:end])
    B = Float64.(Matrix(B_df[!, 2:end]) .> 0)
    S = Float64.(Matrix(S_df[!, 2:end]) .> 0)

    # Presence of each plant and bird in each layer, plants first
    in_B = vcat(vec(sum(B, dims=2)) .> 0, vec(sum(B, dims=1)) .> 0)
    in_S = vcat(vec(sum(S, dims=2)) .> 0, vec(sum(S, dims=1)) .> 0)
    active = in_B .| in_S
    category = [b && s ? "Multilayer" : (b ? "Pollination" : "Seed-dispersal")
                for (b, s) in zip(in_B, in_S)]

    np = length(plants)
    keep_plants = active[1:np]
    keep_birds = active[np+1:end]

    return (network=network,
            species=vcat(plants, birds)[active],
            category=category[active],
            B=B[keep_plants, keep_birds],
            S=S[keep_plants, keep_birds])
end

#--------------------------------------------------------------#
#-Hurricane damage and the empirical regimes-------------------#
#--------------------------------------------------------------#

# Damage fraction removed from every species by a hurricane of each
# Saffir-Simpson category. Damage rises steeply from category 3 to
# categories 4 and 5.
const DAMAGE = Dict("Cat_1" => 0.1, "Cat_2" => 0.2, "Cat_3" => 0.3,
                    "Cat_4" => 0.6, "Cat_5" => 0.8)

# Read the empirical hurricane regime of each island: hurricane frequency
# (lambda) and PMF of categories, with the damage fraction of each category.
# Returns a Dict from island code to a NamedTuple with network, lambda and
# pmf. pmf maps each category to its probability (prob) and damage fraction
# (damage).
function empirical_regimes()
    lambdas = CSV.read(srcdir("empirical_lambdas.csv"), DataFrame)
    pmfs = CSV.read(srcdir("empirical_pmfs.csv"), DataFrame)

    regimes = Dict{String, NamedTuple}()
    for row in eachrow(lambdas)
        island = String(row.network)
        island_pmf = pmfs[pmfs.network .== island, :]
        pmf = Dict(String(r.Category) => (prob=Float64(r.prob), damage=DAMAGE[String(r.Category)])
                   for r in eachrow(island_pmf))
        regimes[island] = (network=island, lambda=Float64(row.lambda), pmf=pmf)
    end
    return regimes
end

#--------------------------------------------------------------#
#-Interaction matrix of the population dynamics----------------#
#--------------------------------------------------------------#

# Add the intra-guild competition of one layer to Y. Every pair of species
# of the guild present in the layer competes with strength c divided by the
# number of other species of the guild in the layer. The denominator is 1
# for a guild with a single species. present marks the species of the guild
# in the layer and offset shifts them to their rows in Y (0 for plants, the
# number of plants for birds).
function add_competition!(; Y::Matrix{Float64}, present::AbstractVector{Bool}, offset::Int64, c::Float64)
    species = findall(present)
    n_others = max(length(species) - 1, 1)
    @inbounds for j in species, i in species
        if i != j
            Y[offset+i, offset+j] += c / n_others
        end
    end
    return nothing
end

# Add the mutualistic benefits of one layer to Y. A plant and a bird that
# interact in A benefit each other with the base strength of the layer
# divided by the degree of the receiving species raised to delta (generalist
# trade-off). With delta = 0 every partner brings the full base strength.
# The higher delta, the less a species gains from each additional partner.
# With delta = 1 the benefits of a species sum to the same total whatever
# its degree. k_plant and k_bird are the degrees used in the trade-off, in
# the order of the rows and columns of A.
function add_mutualism!(; Y::Matrix{Float64}, A::Matrix{Float64}, strength::Float64, k_plant::Vector{Float64}, k_bird::Vector{Float64}, delta::Float64)
    Np = size(A, 1)
    @inbounds for i in 1:Np, j in 1:size(A, 2)
        if A[i, j] > 0.0
            Y[i, Np+j] += strength * A[i, j] / (k_plant[i])^delta
            Y[Np+j, i] += strength * A[i, j] / (k_bird[j])^delta
        end
    end
    return nothing
end

# Build the interaction matrix Y of one multilayer network. Rows and columns
# list the plants first and then the birds, and Y[i, j] is the per-capita
# effect of species j on species i. Each layer adds the competition within
# each guild and the benefits between plants and birds. A pair present in
# both layers competes in both, and a dual-role bird benefits from both.
# The generalist trade-off divides the benefits of a species by its degree
# raised to delta, and tradeoff sets the degree. With "layer" it is the
# number of partners in the layer of the interaction. With "multilayer" it
# is the number of interactions over the two layers. With delta = 1 the
# benefits of a species sum to the base strength in each layer where it
# interacts ("layer") or over the two layers ("multilayer", for m equal to
# d).
# c is negative (-0.1 in the simulations), because the dynamics separate
# competition from mutualism by the sign of the entries of Y. The
# simulations use delta = 0.5 and tradeoff = "layer". Fills the
# preallocated Y in place.
function yij_multi!(; B::Matrix{Float64}, S::Matrix{Float64}, Y::Matrix{Float64}, m::Float64, c::Float64, d::Float64, delta::Float64, tradeoff::String)
    if !(tradeoff in ("layer", "multilayer"))
        error("tradeoff must be \"layer\" or \"multilayer\".")
    end

    Np = size(B, 1)
    fill!(Y, 0.0)

    # Degree of each plant and bird over the two layers
    k_plant_multi = vec(sum(B, dims=2)) .+ vec(sum(S, dims=2))
    k_bird_multi = vec(sum(B, dims=1)) .+ vec(sum(S, dims=1))

    for (A, strength) in ((B, m), (S, d))
        plant_in_layer = vec(sum(A .> 0.0, dims=2) .> 0.0)
        bird_in_layer = vec(sum(A .> 0.0, dims=1) .> 0.0)
        add_competition!(Y=Y, present=plant_in_layer, offset=0, c=c)
        add_competition!(Y=Y, present=bird_in_layer, offset=Np, c=c)

        # Degree used in the trade-off: within the layer or over both layers
        k_plant = tradeoff == "layer" ? vec(sum(A, dims=2)) : k_plant_multi
        k_bird = tradeoff == "layer" ? vec(sum(A, dims=1)) : k_bird_multi
        add_mutualism!(Y=Y, A=A, strength=strength, k_plant=k_plant, k_bird=k_bird, delta=delta)
    end

    return nothing
end
