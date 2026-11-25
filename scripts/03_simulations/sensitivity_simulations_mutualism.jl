#--------------------------------------------------------------#
#-Sensitivity of the simulations to the mutualistic strength---#
#--------------------------------------------------------------#

# Repeat the main simulations with weaker and stronger mutualistic benefits.
# The base strength of pollination and of seed dispersal (m and d) is 0.5 in
# the main simulations. Here both take 0.25, 0.35, 0.75 and 1.0. For each
# strength, run the frequency scenario (880 combinations of island and
# lambda), the intensity scenario (231 combinations of island and shift)
# and the null networks (4,400). Each combination runs 1,000 replicates of
# 300 years. The four strengths give 22,044 combinations of parameters.
# Edit experiments, null_models and strengths to run fewer.
#
# Inputs:  src/multi_pol/B_<island>.csv
#          src/multi_sd/S_<island>.csv
#          src/null_<model>/network_manifest.csv and the matrices it lists
#          src/empirical_lambdas.csv
#          src/empirical_pmfs.csv
# Outputs: data/sens_sims_lambda_gamma<strength>/*.csv
#          data/sens_sims_pmf_gamma<strength>/*.csv
#          data/sens_sims_null_gamma<strength>/*.csv
#          (files named as in the main simulations)

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
    #-1. Experiments, strengths and parameters---------------------#
    #--------------------------------------------------------------#

    # Experiments to run: "frequency" and "intensity" (the two hurricane
    # scenarios of the empirical networks) and "null" (null networks under
    # the empirical regimes)
    experiments = ["frequency", "intensity", "null"]

    # Null models of the null experiment. Null model X is in src/null_X/
    null_models = ["I", "I-A", "II", "II-B"]

    # Base mutualistic strengths, the same for pollination and seed
    # dispersal. With the 0.5 of the main simulations they form a gradient
    # of five strengths, from half to twice the main one. The saturating
    # functional response keeps the dynamics bounded at any strength, and in
    # this range every network reaches a steady state before the hurricanes.
    # With much weaker benefits, species with low growth rates decline
    # slowly instead of settling, and the extinction threshold (1% of the
    # steady-state abundance) loses its reference.
    strengths = [0.25, 0.35, 0.75, 1.0]

    # Parameters shared by all scenarios, as in the main simulations
    # (generalist trade-off with delta = 0.5 within each layer)
    params = Dict(
        :c => -0.1,
        :ht => 0.5,
        :delta => 0.5,
        :tradeoff => "layer",
        :t_years => 300,
        :nsim => 1000
    )

    #--------------------------------------------------------------#
    #-2. Empirical networks and regimes----------------------------#
    #--------------------------------------------------------------#

    # One network per island, named after the island code of its files
    B_paths = readdir(srcdir("multi_pol"); join=true)
    S_paths = readdir(srcdir("multi_sd"); join=true)
    networks = [replace(basename(path), r"^B_|\.csv$" => "") for path in B_paths]

    net_data = [load_network(B_df=CSV.read(B_paths[i], DataFrame),
                             S_df=CSV.read(S_paths[i], DataFrame),
                             network=networks[i]) for i in eachindex(networks)]

    # Empirical lambda and PMF of each island, by island code
    regimes = empirical_regimes()

    # Scenarios of all experiments, one dictionary per simulation file
    c_list = Dict{Symbol, Any}[]

    #--------------------------------------------------------------#
    #-3. Frequency scenario: island, lambda and strength-----------#
    #--------------------------------------------------------------#

    if "frequency" in experiments

        # Values of lambda, in hurricanes per year
        lambda_val = collect(0.01:0.01:0.80)

        # Each island with every value of lambda and its own empirical PMF,
        # for each strength
        c_lambda = vec([merge(params, Dict(
            :net => net,
            :m => strength,
            :d => strength,
            :lambda => lambda,
            :intensity_scenario => regimes[net.network].pmf,
            :out_dir => "sens_sims_lambda_gamma$(strength)"
        )) for net in net_data, lambda in lambda_val, strength in strengths])

        append!(c_list, c_lambda)
    end

    #--------------------------------------------------------------#
    #-4. Intensity scenario: island, shift and strength------------#
    #--------------------------------------------------------------#

    # Move probability mass from categories 1 to 3 to categories 4 and 5.
    # shift is the share of probability mass moved, at most the mass of
    # categories 1 to 3 (a shift of 1 moves all of it). Within each group the
    # categories keep their empirical proportions.
    # Returns a Dict like pmf, with the shifted probabilities.
    function shift_pmf(; pmf::Dict, shift::Float64)
        low = ("Cat_1", "Cat_2", "Cat_3")
        p_low = sum(pmf[cat].prob for cat in low if haskey(pmf, cat); init=0.0)
        p_high = sum(pmf[cat].prob for cat in ("Cat_4", "Cat_5") if haskey(pmf, cat); init=0.0)
        moved = min(shift, p_low)

        # Categories 1 to 3 lose the mass moved and categories 4 and 5 gain
        # it, each in proportion to its probability
        scale_low = p_low > 0.0 ? (p_low - moved) / p_low : 0.0
        scale_high = p_high > 0.0 ? (p_high + moved) / p_high : 0.0

        return Dict(cat => (prob=v.prob * (cat in low ? scale_low : scale_high), damage=v.damage)
                    for (cat, v) in pmf)
    end

    if "intensity" in experiments

        # Values of the shift. Each step moves another 5% of the probability
        # mass to categories 4 and 5. A shift of 0 is the empirical regime.
        shifts = collect(0.0:0.05:1.0)

        # Each island with its own empirical lambda and its shifted PMF, for
        # each strength
        c_pmf = vec([merge(params, Dict(
            :net => net,
            :m => strength,
            :d => strength,
            :lambda => regimes[net.network].lambda,
            :intensity_scenario => shift_pmf(pmf=regimes[net.network].pmf, shift=shift),
            :shift => shift,
            :out_dir => "sens_sims_pmf_gamma$(strength)"
        )) for net in net_data, shift in shifts, strength in strengths])

        append!(c_list, c_pmf)
    end

    #--------------------------------------------------------------#
    #-5. Null networks: null network and strength------------------#
    #--------------------------------------------------------------#

    if "null" in experiments

        # Stack the manifests of the null models: one row per null network,
        # with its model, island, replicate and matrix files
        null_manifest = vcat([CSV.read(srcdir("null_$(model)", "network_manifest.csv"), DataFrame)
                              for model in null_models]...; cols=:union)

        # Every null network with the code of its island. Network name:
        # <island>_null_<model>_<id>, with the id padded to four digits, as
        # in data/sims/.
        null_data = [(island=String(row.network),
                      net=load_network(B_df=CSV.read(srcdir("null_$(row.model)", row.pol_file), DataFrame),
                                       S_df=CSV.read(srcdir("null_$(row.model)", row.sd_file), DataFrame),
                                       network="$(row.network)_null_$(row.model)_$(lpad(string(row.null_id), 4, '0'))"))
                     for row in eachrow(null_manifest)]

        # Each null network with the empirical lambda and PMF of its island,
        # for each strength
        c_null = vec([merge(params, Dict(
            :net => null_net.net,
            :m => strength,
            :d => strength,
            :lambda => regimes[null_net.island].lambda,
            :intensity_scenario => regimes[null_net.island].pmf,
            :out_dir => "sens_sims_null_gamma$(strength)"
        )) for null_net in null_data, strength in strengths])

        append!(c_list, c_null)
    end

    println("Scenarios to run: $(length(c_list))")

end

#--------------------------------------------------------------#
#-6. Run the simulations---------------------------------------#
#--------------------------------------------------------------#

pmap(net_multisim_typeII, c_list)
