#--------------------------------------------------------------#
#-Empirical networks under increased hurricane frequency-------#
#--------------------------------------------------------------#

# Simulate the 11 empirical multilayer networks under increased hurricane
# frequency. Each island keeps its empirical PMF of hurricane categories and
# lambda goes from 0.01 to 0.80 hurricanes per year in steps of 0.01. The
# 880 combinations of island and lambda each run 1,000 replicates of 300
# years.
#
# Inputs:  src/multi_pol/B_<island>.csv
#          src/multi_sd/S_<island>.csv
#          src/empirical_lambdas.csv
#          src/empirical_pmfs.csv
# Outputs: data/main_sims_lambda/
#          <island>_multilayer_lambda<lambda>_pmf<d>_tyears300.csv
#          (d is the expected damage fraction of the empirical PMF)

using Distributed

# Worker processes. pmap gives one scenario at a time to each worker.
addprocs(12)

@everywhere using DrWatson
@everywhere @quickactivate "Goncalves_et_al_shifting_hurricanes"

@everywhere begin

    using DifferentialEquations
    using Distributions
    using LinearAlgebra
    using DrWatson
    using CSV
    using DataFrames
    using Random

    BLAS.set_num_threads(1)

    # Functions of the population dynamics, networks and hurricane regimes
    include(srcdir("main_functions.jl"))
    include(srcdir("aux_functions.jl"))

    #--------------------------------------------------------------#
    #-1. Empirical networks----------------------------------------#
    #--------------------------------------------------------------#

    # One network per island, named after the island code of its files
    B_paths = readdir(srcdir("multi_pol"); join=true)
    S_paths = readdir(srcdir("multi_sd"); join=true)
    networks = [replace(basename(path), r"^B_|\.csv$" => "") for path in B_paths]

    net_data = [load_network(B_df=CSV.read(B_paths[i], DataFrame),
                             S_df=CSV.read(S_paths[i], DataFrame),
                             network=networks[i]) for i in eachindex(networks)]

    #--------------------------------------------------------------#
    #-2. Empirical regimes and values of lambda--------------------#
    #--------------------------------------------------------------#

    # Empirical regime of each island. Only the PMF is used here.
    regimes = collect(values(empirical_regimes()))

    # Values of lambda, in hurricanes per year
    lambda_val = collect(0.01:0.01:0.80)

    #--------------------------------------------------------------#
    #-3. One scenario per island and value of lambda---------------#
    #--------------------------------------------------------------#

    c1 = Dict(
        :net => net_data,
        :regime => regimes,
        :lambda => lambda_val,
        :m => 0.5,
        :d => 0.5,
        :c => -0.1,
        :ht => 0.5,
        :delta => 0.5,
        :tradeoff => "layer",
        :t_years => 300,
        :nsim => 1000,
        :out_dir => "main_sims_lambda"
    )

    # All combinations of network, regime and lambda. Keep those with the
    # regime of the island of the network, and pass its PMF to the wrapper.
    c_list = filter(d -> d[:net].network == d[:regime].network, dict_list(c1))
    for d in c_list
        d[:intensity_scenario] = d[:regime].pmf
        delete!(d, :regime)
    end

    println("Scenarios to run: $(length(c_list))")

end

#--------------------------------------------------------------#
#-4. Run the simulations---------------------------------------#
#--------------------------------------------------------------#

pmap(net_multisim_typeII, c_list)
