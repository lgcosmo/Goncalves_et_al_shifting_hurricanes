#--------------------------------------------------------------#
#-Null networks under the empirical hurricane regime-----------#
#--------------------------------------------------------------#

# Simulate the population dynamics of every null network of the decoupling
# null models I, I-A, II and II-B under the empirical hurricane regime of
# its island. Each null network keeps its island's empirical hurricane
# frequency (lambda) and PMF of categories, and the parameters of the
# empirical simulations. Each null model has 100 networks per island on 11
# islands, so the four models schedule 4,400 tasks of 1,000 replicates of
# 300 years. The extinctions of the null networks are compared with those of
# the empirical networks in Extended Data Fig. 1 (process_null_simulations.R
# and null_model_comparison.R). Edit null_models below to run fewer models.
#
# Inputs:  src/null_<model>/network_manifest.csv and the matrices it lists,
#          built by scripts/01_network_preparation/build_null_decoupling.jl
#          src/empirical_lambdas.csv and src/empirical_pmfs.csv, written by
#          empirical_lambda_pmf.R
# Outputs: data/sims/<island>_null_<model>_<id>_multilayer_lambda<lambda>_pmf<d>_tyears300.csv,
#          one file per null network, read by process_null_simulations.R.
#          The network name identifies the island, the null model and the
#          null network, so the four models never overwrite each other.
# Helpers: src/main_functions.jl (population dynamics, hurricanes and the
#          simulation wrapper) and src/aux_functions.jl (networks, regimes
#          and the interaction matrix)

using Distributed

# Worker processes; pmap gives each one scenario at a time
addprocs(12)

@everywhere using DrWatson
@everywhere @quickactivate "Goncalves_et_al_shifting_hurricanes"

@everywhere begin

    using DifferentialEquations
    using Distributions
    using LinearAlgebra
    using CSV
    using DataFrames
    using Random
    using Statistics

    BLAS.set_num_threads(1)

    # Population dynamics and the simulation wrapper, then the networks, the
    # regimes and the interaction matrix
    include(srcdir("main_functions.jl"))
    include(srcdir("aux_functions.jl"))

    #--------------------------------------------------------------#
    #-1. Null networks of the four null models--------------------#
    #--------------------------------------------------------------#

    # Null models to simulate. Each null model X lives in src/null_X/
    null_models = ["I", "I-A", "II", "II-B"]

    # Stack the manifests of all null models, so that every null network is
    # one row with its model, island, replicate and matrix files
    null_manifest = vcat([CSV.read(srcdir("null_$(model)", "network_manifest.csv"), DataFrame)
                          for model in null_models]...; cols=:union)

    # The network name carries the island, the null model and the null
    # network, which process_null_simulations.R parses from the file name
    net_data = [load_network(B_df=CSV.read(srcdir("null_$(row.model)", row.pol_file), DataFrame),
                             S_df=CSV.read(srcdir("null_$(row.model)", row.sd_file), DataFrame),
                             network="$(row.network)_null_$(row.model)_$(lpad(string(row.null_id), 4, '0'))")
                for row in eachrow(null_manifest)]

    #--------------------------------------------------------------#
    #-2. Empirical hurricane regime of each island-----------------#
    #--------------------------------------------------------------#

    # Frequency and PMF of each island, so that every null network meets the
    # regime of its own island
    regimes = collect(values(empirical_regimes()))

    #--------------------------------------------------------------#
    #-3. One scenario per null network----------------------------#
    #--------------------------------------------------------------#

    c1 = Dict(
        :net => net_data,
        :regime => regimes,
        :m => 0.5,
        :d => 0.5,
        :c => -0.1,
        :ht => 0.5,
        :t_years => 300,
        :nsim => 1000
    )

    # Every combination of null network and island regime, then only the
    # regime of the null network's own island, whose frequency and PMF the
    # wrapper expects
    c_list = filter(d -> startswith(d[:net].network, d[:regime].network * "_null_"), dict_list(c1))
    for d in c_list
        d[:lambda] = d[:regime].lambda
        d[:intensity_scenario] = d[:regime].pmf
        delete!(d, :regime)
    end

end

#--------------------------------------------------------------#
#-4. Run the simulations---------------------------------------#
#--------------------------------------------------------------#

println("Null networks to simulate: $(length(c_list))")

pmap(net_multisim_typeII, c_list)
