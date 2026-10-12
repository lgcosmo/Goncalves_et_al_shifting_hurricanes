#--------------------------------------------------------------#
#-Reading and writing the multilayer networks of the islands---#
#--------------------------------------------------------------#

# Helpers shared by null models I to II-B, the decoupling null models of
# Extended Data Fig. 1. They read the empirical pollination and
# seed-dispersal matrices of each island and the list of hummingbird species,
# derive a reproducible seed for each null network, and write null matrices
# in the same format as the empirical ones. Plants and birds keep their
# labels and order, and positive observation counts become binary
# interactions, as in the loader of the simulations (load_network in
# src/aux_functions.jl).
#
# Inputs:  src/multi_pol/B_<island>.csv (plants x birds, flower visitation)
#          src/multi_sd/S_<island>.csv  (plants x birds, fruit consumption)
#          src/hummingbirds.csv (one hummingbird species per row)
# Outputs: the null matrices of src/null_<model>/, written through
#          write_null_matrix by build_null_ensemble
# Used by: src/null_decoupling_functions.jl, included after this file by
#          scripts/01_network_preparation/build_null_decoupling.jl

using CSV
using DataFrames
using SHA

#--------------------------------------------------------------#
#-Reading the empirical networks-------------------------------#
#--------------------------------------------------------------#

# Normalize a bird name only to match it against the hummingbird list. Matrix
# labels are never replaced by the normalized names.
# Arguments: name is a bird common name. Returns the normalized String.
function normalize_bird_name(; name::AbstractString)
    cleaned = replace(String(name), r"[‐‑‒–—−]" => "-")
    return lowercase(strip(replace(cleaned, r"\s+" => " ")))
end

# Read the hummingbird species, the birds that never receive a fruit
# interaction in the null networks.
# Arguments: path is the CSV file with one species per row in the column
# species. Returns the normalized names as a Set.
function read_hummingbirds(; path::String)
    names = String.(CSV.read(path, DataFrame).species)
    isempty(names) && error("No hummingbird names found in $path")
    return Set(normalize_bird_name(name=name) for name in names)
end

# Mark the hummingbirds among the birds of a network.
# Arguments: birds holds the bird labels of the network; hummingbirds is the
# Set returned by read_hummingbirds. Returns a Bool vector, one per bird.
function hummingbird_mask(; birds::Vector{String}, hummingbirds::Set{String})
    return Bool[normalize_bird_name(name=bird) in hummingbirds for bird in birds]
end

# Read one incidence matrix without changing its labels or order. The first
# column, named plant, holds the plant labels; the headers hold the birds.
# Arguments: path is the CSV file. Returns the labels and the binary matrix.
function read_incidence_matrix(; path::String)
    df = CSV.read(path, DataFrame; normalizenames=false)
    names(df)[1] == "plant" || error("First column must be named plant: $path")
    plants = String.(df[!, 1])
    birds = String.(names(df)[2:end])
    length(unique(plants)) == length(plants) || error("Duplicate plant names: $path")
    length(unique(birds)) == length(birds) || error("Duplicate bird names: $path")

    # Observation counts become binary interactions
    counts = Matrix(df[:, 2:end])
    all(x -> x isa Real && isfinite(x) && x >= 0, counts) ||
        error("Matrix entries must be finite, nonnegative numbers: $path")
    return (plants=plants, birds=birds, matrix=Int64.(counts .> 0))
end

# Read the pollination (B) and seed-dispersal (S) matrices of every island.
# Both layers of an island must list the same plants and birds, in the same
# order. Arguments: source_dir is src/. Returns one NamedTuple per island with
# its name, labels, binary layers and the SHA-256 of its two input files.
function read_empirical_networks(; source_dir::String)
    pol_dir = joinpath(source_dir, "multi_pol")
    sd_dir = joinpath(source_dir, "multi_sd")
    islands = sort([name[3:end-4] for name in readdir(pol_dir)
                    if startswith(name, "B_") && endswith(name, ".csv")])
    isempty(islands) && error("No empirical pollination matrices found.")

    networks = NamedTuple[]
    for island in islands
        b_path = joinpath(pol_dir, "B_$(island).csv")
        s_path = joinpath(sd_dir, "S_$(island).csv")
        pol = read_incidence_matrix(path=b_path)
        sd = read_incidence_matrix(path=s_path)
        pol.plants == sd.plants && pol.birds == sd.birds ||
            error("Species labels or order differ between layers for $island")
        push!(networks, (network=island, plants=pol.plants, birds=pol.birds,
                         B=pol.matrix, S=sd.matrix,
                         pol_sha256=bytes2hex(sha256(read(b_path))),
                         sd_sha256=bytes2hex(sha256(read(s_path)))))
    end
    return networks
end

#--------------------------------------------------------------#
#-Seeds and writing the null networks--------------------------#
#--------------------------------------------------------------#

# Derive an independent, reproducible seed for one null network from the base
# seed of its null model, its island and its replicate number. SHA-256 avoids
# Julia's session-dependent hashing.
# Arguments identify the null network. Returns a UInt64 seed for Xoshiro.
function null_network_seed(; base_seed::Int64, network::String, null_id::Int64)
    digest = bytes2hex(sha256("$base_seed|$network|$null_id"))
    return parse(UInt64, digest[1:16]; base=16)
end

# Write one binary matrix in the empirical format: the plant column first,
# then one column per bird, with the original labels and order.
# Arguments: matrix is the layer; plants and birds are its labels; path is the
# destination, which must not exist yet. Returns nothing.
function write_null_matrix(; matrix::Matrix{Int64}, plants::Vector{String},
                             birds::Vector{String}, path::String)
    ispath(path) && error("Refusing to overwrite existing file: $path")
    df = DataFrame(matrix, birds)
    insertcols!(df, 1, :plant => plants)
    CSV.write(path, df; quotestrings=true)
    return nothing
end
