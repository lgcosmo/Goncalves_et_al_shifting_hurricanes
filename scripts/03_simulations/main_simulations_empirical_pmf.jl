#--------------------------------------------------------------#
#-Empirical networks under increased hurricane intensity-------#
#--------------------------------------------------------------#

# Simulate the 11 empirical multilayer networks under increased hurricane
# intensity. Each island keeps its empirical lambda, and a shift s moves a
# share of the probability mass of its PMF from categories 1 to 3 to
# categories 4 and 5. The shift goes from 0 to 1 in steps of 0.05 (21
# values). The 231 combinations of island and shift each run 1,000
# replicates of 300 years. A shift of 0 is the empirical regime of the
# island, the baseline of the null-model comparison.
#
# Inputs:  src/multi_pol/B_<island>.csv
#          src/multi_sd/S_<island>.csv
#          src/empirical_lambdas.csv
#          src/empirical_pmfs.csv
# Outputs: data/main_sims_pmf/
#          <island>_multilayer_lambda<lambda>_pmf<s>_tyears300.csv

using Distributed

# Worker processes. pmap gives one scenario at a time to each worker.
addprocs(16)

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
    #-2. Empirical regime of each island---------------------------#
    #--------------------------------------------------------------#

    # Empirical lambda and PMF of each island, by island code
    regimes = empirical_regimes()

    #--------------------------------------------------------------#
    #-3. Shift of probability mass toward categories 4 and 5-------#
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

    #--------------------------------------------------------------#
    #-4. One scenario per island and shift-------------------------#
    #--------------------------------------------------------------#

    # Values of the shift. Each step moves another 5% of the probability
    # mass to categories 4 and 5. A shift of 0 is the empirical regime.
    shifts = collect(0.0:0.05:1.0)

    # Each island with its own empirical lambda and its shifted PMF
    c_list = vec([Dict(
        :net => net,
        :lambda => regimes[net.network].lambda,
        :intensity_scenario => shift_pmf(pmf=regimes[net.network].pmf, shift=shift),
        :shift => shift,
        :m => 0.5,
        :d => 0.5,
        :c => -0.1,
        :ht => 0.5,
        :delta => 0.5,
        :tradeoff => "layer",
        :t_years => 300,
        :nsim => 1000,
        :out_dir => "main_sims_pmf"
    ) for net in net_data, shift in shifts])

    println("Scenarios to run: $(length(c_list))")

end

#--------------------------------------------------------------#
#-5. Run the simulations---------------------------------------#
#--------------------------------------------------------------#

pmap(net_multisim_typeII, c_list)
