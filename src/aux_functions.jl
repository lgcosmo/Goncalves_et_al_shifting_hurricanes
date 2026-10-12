#--------------------------------------------------------------#
#-Networks, hurricane regimes and the interaction matrix-------#
#--------------------------------------------------------------#

# Functions that prepare the inputs of the simulations. load_network reads
# the two incidence matrices of one island or null network and lists its
# species. DAMAGE holds the damage fraction of each hurricane category and
# empirical_regimes reads the empirical hurricane regime of every island.
# yij_multi! builds the interaction matrix of the population dynamics model
# of the Methods.
#
# Inputs:  the incidence matrices of a network, read as DataFrames by the
#          simulation scripts from src/multi_pol, src/multi_sd or
#          src/null_<model>
#          src/empirical_lambdas.csv and src/empirical_pmfs.csv
#          (empirical_lambda_pmf.R)
# Outputs: none on disk. The network and the regime feed net_multisim_typeII
#          and the interaction matrix feeds eco_dynamics_typeII!, both in
#          src/main_functions.jl.
# Used by: the three simulation scripts in scripts/03_simulations, which
#          load CSV, DataFrames and DrWatson and include this file after
#          main_functions.jl

using CSV
using DataFrames
using LinearAlgebra

#--------------------------------------------------------------#
#-Species and binary layers of one network---------------------#
#--------------------------------------------------------------#

# Read one multilayer network from its two incidence matrices. Every positive
# entry becomes an interaction, and only the species with at least one
# interaction in either layer are kept, plants first and then birds. A
# species present in both layers is "Multilayer" (for birds, the dual-role
# birds), and the others are "Pollination" or "Seed-dispersal". This category
# is the species type of the simulation files.
# Arguments: B_df and S_df are the pollination and seed-dispersal matrices
# read as DataFrames, with the plants in the first column and one column per
# bird; network is the name of the network, the island code or the null
# network id. Returns a NamedTuple with network, species (their names),
# category, and the binary matrices B and S restricted to the species kept.
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

# Damage fraction that a hurricane of each Saffir-Simpson category removes
# from every species, as in the Methods. Losses rise steeply from category 3
# to categories 4 and 5, the catastrophic storms of the record.
const DAMAGE = Dict("Cat_1" => 0.1, "Cat_2" => 0.2, "Cat_3" => 0.3,
                    "Cat_4" => 0.6, "Cat_5" => 0.8)

# Read the empirical hurricane regime of every island, estimated from the
# historical record by empirical_lambda_pmf.R. The regime combines the
# hurricane frequency lambda and the PMF of categories, with the damage
# fraction of each category attached.
# Arguments: none. Returns a Dict from island code to a NamedTuple with
# network (the island code), lambda, and pmf, which maps each category to a
# NamedTuple with its probability (prob) and damage fraction (damage).
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

# Add the intra-guild competition of one layer to the interaction matrix.
# Every pair of species of the guild present in the layer competes with
# strength c divided by the number of other species of the guild in that
# layer. The denominator stays at 1 when the guild has a single species.
# Arguments: Y is the interaction matrix; present marks the species of the
# guild present in the layer, in the order of their rows in Y after offset
# (0 for plants, the number of plants for birds); c is the competition
# strength. Fills Y in place and returns nothing.
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

# Add the mutualistic benefits of one layer to the interaction matrix. A
# plant and a bird that interact benefit each other with the base strength
# of the layer divided by the square root of the number of partners of the
# receiving species in that layer, the generalist trade-off with exponent
# delta = 0.5.
# Arguments: Y is the interaction matrix; A is the binary incidence matrix
# of the layer (plants in rows, birds in columns); strength is the base
# mutualistic strength of the layer. Fills Y in place and returns nothing.
function add_mutualism!(; Y::Matrix{Float64}, A::Matrix{Float64}, strength::Float64)
    Np = size(A, 1)
    k_plant = vec(sum(A, dims=2))
    k_bird = vec(sum(A, dims=1))
    @inbounds for i in 1:Np, j in 1:size(A, 2)
        if A[i, j] > 0.0
            Y[i, Np+j] += strength * A[i, j] / (k_plant[i])^0.5
            Y[Np+j, i] += strength * A[i, j] / (k_bird[j])^0.5
        end
    end
    return nothing
end

# Build the interaction matrix Y of the population dynamics model of the
# Methods for one multilayer network. Rows and columns list the plants first
# and then the birds, and Y[i, j] is the per-capita effect of species j on
# species i. Each layer adds the competition within each guild and the
# mutualistic benefits between plants and birds, so a pair present in both
# layers competes in both and a dual-role bird benefits from both.
# Arguments: B and S are the binary incidence matrices of the pollination and
# seed-dispersal layers (plants in rows, birds in columns); Y is the
# preallocated square matrix with one row per plant and bird; m and d are the
# base mutualistic strengths of pollination and seed dispersal; c is the
# intra-guild competition strength, negative (-0.1 in the simulations), so
# that the dynamics can separate the competitive from the mutualistic entries
# by their sign. Fills Y in place and returns nothing.
function yij_multi!(; B::Matrix{Float64}, S::Matrix{Float64}, Y::Matrix{Float64}, m::Float64, c::Float64, d::Float64)
    Np = size(B, 1)
    fill!(Y, 0.0)

    for (A, strength) in ((B, m), (S, d))
        plant_in_layer = vec(sum(A .> 0.0, dims=2) .> 0.0)
        bird_in_layer = vec(sum(A .> 0.0, dims=1) .> 0.0)
        add_competition!(Y=Y, present=plant_in_layer, offset=0, c=c)
        add_competition!(Y=Y, present=bird_in_layer, offset=Np, c=c)
        add_mutualism!(Y=Y, A=A, strength=strength)
    end

    return nothing
end
