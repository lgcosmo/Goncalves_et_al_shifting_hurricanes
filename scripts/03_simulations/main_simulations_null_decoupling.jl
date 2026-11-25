#--------------------------------------------------------------#
#-Null networks under the empirical hurricane regime-----------#
#--------------------------------------------------------------#

# Simulate the null networks of the decoupling null models I, I-A, II and
# II-B under the empirical hurricane regime of their island (lambda and PMF
# of categories), with the parameters of the empirical simulations. Each
# null model has 100 networks per island on 11 islands. The four models give
# 4,400 tasks of 1,000 replicates of 300 years. Edit null_models to run
# fewer models.
#
# Inputs:  src/null_<model>/network_manifest.csv and the matrices it lists
#          src/empirical_lambdas.csv
#          src/empirical_pmfs.csv
# Outputs: data/sims/<network>_multilayer_lambda<lambda>_pmf<d>_tyears300.csv
#          (<network> is <island>_null_<model>_<id> and d the expected
#          damage fraction of the empirical PMF)

using Distributed

# Worker processes. pmap gives one scenario at a time to each worker.
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

    # Functions of the population dynamics, networks and hurricane regimes
    include(srcdir("main_functions.jl"))
    include(srcdir("aux_functions.jl"))

    #--------------------------------------------------------------#
    #-1. Null networks of the four null models---------------------#
    #--------------------------------------------------------------#

    # Null models to simulate. Null model X is in src/null_X/
    null_models = ["I", "I-A", "II", "II-B"]

    # Stack the manifests of the null models: one row per null network, with
    # its model, island, replicate and matrix files
    null_manifest = vcat([CSV.read(srcdir("null_$(model)", "network_manifest.csv"), DataFrame)
                          for model in null_models]...; cols=:union)

    # Network name: <island>_null_<model>_<id>, with the id padded to four
    # digits. The four models never overwrite the files of each other.
    net_data = [load_network(B_df=CSV.read(srcdir("null_$(row.model)", row.pol_file), DataFrame),
                             S_df=CSV.read(srcdir("null_$(row.model)", row.sd_file), DataFrame),
                             network="$(row.network)_null_$(row.model)_$(lpad(string(row.null_id), 4, '0'))")
                for row in eachrow(null_manifest)]

    #--------------------------------------------------------------#
    #-2. Empirical hurricane regime of each island-----------------#
    #--------------------------------------------------------------#

    # Empirical lambda and PMF of each island
    regimes = collect(values(empirical_regimes()))

    #--------------------------------------------------------------#
    #-3. One scenario per null network-----------------------------#
    #--------------------------------------------------------------#

    c1 = Dict(
        :net => net_data,
        :regime => regimes,
        :m => 0.5,
        :d => 0.5,
        :c => -0.1,
        :ht => 0.5,
        :delta => 0.5,
        :tradeoff => "layer",
        :t_years => 300,
        :nsim => 1000
    )

    # All combinations of null network and regime. Keep those with the
    # regime of the island of the null network, and pass its lambda and PMF
    # to the wrapper.
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
