#--------------------------------------------------------------#
#-Reading and writing the multilayer networks of the islands---#
#--------------------------------------------------------------#

# Functions shared by null models I, I-A, II and II-B. Read the empirical
# pollination and seed-dispersal matrices of each island and the list of
# hummingbirds, derive a reproducible seed for each null network, and write
# the null matrices in the format of the empirical ones. Plants and birds
# keep their labels and order. Positive observation counts become binary
# interactions, as in load_network (src/aux_functions.jl).
#
# Inputs:  src/multi_pol/B_<island>.csv
#          src/multi_sd/S_<island>.csv
#          src/hummingbirds.csv
# Outputs: src/null_<model>/multi_pol/B_<island>_null_<model>_<id>.csv
#          src/null_<model>/multi_sd/S_<island>_null_<model>_<id>.csv

using CSV
using DataFrames
using SHA

#--------------------------------------------------------------#
#-Reading the empirical networks-------------------------------#
#--------------------------------------------------------------#

# Normalize a bird name (dashes, spaces and case), only to match it against
# the hummingbird list. The labels of the matrices are never replaced.
function normalize_bird_name(; name::AbstractString)
    cleaned = replace(String(name), r"[‐‑‒–—−]" => "-")
    return lowercase(strip(replace(cleaned, r"\s+" => " ")))
end

# Read the hummingbird species (column species of the CSV file). These birds
# never receive a fruit interaction in the null networks. Returns the
# normalized names as a Set.
function read_hummingbirds(; path::String)
    names = String.(CSV.read(path, DataFrame).species)
    isempty(names) && error("No hummingbird names found in $path")
    return Set(normalize_bird_name(name=name) for name in names)
end

# Mark the hummingbirds among the birds of a network (one Bool per bird).
function hummingbird_mask(; birds::Vector{String}, hummingbirds::Set{String})
    return Bool[normalize_bird_name(name=bird) in hummingbirds for bird in birds]
end

# Read one incidence matrix without changing its labels or order. The first
# column (plant) has the plant labels and the headers have the birds.
# Returns the labels and the binary matrix.
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
# Both layers of an island must have the same plants and birds, in the same
# order. Returns one NamedTuple per island with its name, labels, binary
# layers and the SHA-256 of its two input files.
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

# Reproducible seed of one null network, from the base seed of its null
# model, its island and its replicate number. SHA-256 avoids the
# session-dependent hashing of Julia. Returns a UInt64 seed for Xoshiro.
function null_network_seed(; base_seed::Int64, network::String, null_id::Int64)
    digest = bytes2hex(sha256("$base_seed|$network|$null_id"))
    return parse(UInt64, digest[1:16]; base=16)
end

# Write one binary matrix in the empirical format: the plant column first,
# then one column per bird, with the original labels and order. path must
# not exist yet.
function write_null_matrix(; matrix::Matrix{Int64}, plants::Vector{String},
                             birds::Vector{String}, path::String)
    ispath(path) && error("Refusing to overwrite existing file: $path")
    df = DataFrame(matrix, birds)
    insertcols!(df, 1, :plant => plants)
    CSV.write(path, df; quotestrings=true)
    return nothing
end
