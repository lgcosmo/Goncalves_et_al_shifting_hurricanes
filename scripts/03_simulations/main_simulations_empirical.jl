#--------------------------------------------------------------#
#-Empirical networks under increased hurricane frequency-------#
#--------------------------------------------------------------#

# Simulate the population dynamics of the 11 empirical multilayer networks
# under increased hurricane frequency, the first scenario of the Methods.
# Every island keeps its empirical PMF of hurricane categories, and the
# hurricane frequency lambda sweeps from 0.01 to 0.80 hurricanes per year in
# steps of 0.01. Each scenario runs 1,000 replicates of 300 years, so the
# sweep schedules 880 tasks of 1,000 replicates. The extinction curves
# against frequency are Fig. 3c (survival_curves.R and figure_3.R) and give
# the species fragility to hurricane frequency of the pSEM in Fig. 4
# (calculate_auc_extinction.R).
#
# Inputs:  src/multi_pol/B_<island>.csv and src/multi_sd/S_<island>.csv, the
#          incidence matrices of the two layers, built by
#          scripts/01_network_preparation/build_networks_multilayer.R
#          src/empirical_lambdas.csv and src/empirical_pmfs.csv, written by
#          empirical_lambda_pmf.R
# Outputs: data/main_sims_lambda/<island>_multilayer_lambda<lambda>_pmf<d>_tyears300.csv,
#          one file per island and frequency, where d is the expected damage
#          fraction of the island's empirical PMF, read by survival_curves.R
# Helpers: src/main_functions.jl (population dynamics, hurricanes and the
#          simulation wrapper) and src/aux_functions.jl (networks, regimes
#          and the interaction matrix). main_simulations_empirical_pmf.jl
#          runs the companion sweep over hurricane intensity.

using Distributed

# Worker processes; pmap gives each one scenario at a time
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

    # Population dynamics and the simulation wrapper, then the networks, the
    # regimes and the interaction matrix
    include(srcdir("main_functions.jl"))
    include(srcdir("aux_functions.jl"))

    #--------------------------------------------------------------#
    #-1. Empirical networks----------------------------------------#
    #--------------------------------------------------------------#

    # One network per island, named by the island code of its matrix files
    B_paths = readdir(srcdir("multi_pol"); join=true)
    S_paths = readdir(srcdir("multi_sd"); join=true)
    networks = [replace(basename(path), r"^B_|\.csv$" => "") for path in B_paths]

    net_data = [load_network(B_df=CSV.read(B_paths[i], DataFrame),
                             S_df=CSV.read(S_paths[i], DataFrame),
                             network=networks[i]) for i in eachindex(networks)]

    #--------------------------------------------------------------#
    #-2. Empirical regimes and the frequency sweep-----------------#
    #--------------------------------------------------------------#

    # The empirical regime of each island, of which only the PMF is kept
    regimes = collect(values(empirical_regimes()))

    # Hurricane frequencies of the sweep, in hurricanes per year
    lambda_val = collect(0.01:0.01:0.80)

    #--------------------------------------------------------------#
    #-3. One scenario per island and frequency---------------------#
    #--------------------------------------------------------------#

    c1 = Dict(
        :net => net_data,
        :regime => regimes,
        :lambda => lambda_val,
        :m => 0.5,
        :d => 0.5,
        :c => -0.1,
        :ht => 0.5,
        :t_years => 300,
        :nsim => 1000,
        :out_dir => "main_sims_lambda"
    )

    # Every combination of network, regime and frequency, then only the
    # regime of the network's own island, whose PMF the wrapper expects
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
