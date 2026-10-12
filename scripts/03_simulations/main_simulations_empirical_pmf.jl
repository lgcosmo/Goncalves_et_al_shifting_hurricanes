#--------------------------------------------------------------#
#-Empirical networks under increased hurricane intensity-------#
#--------------------------------------------------------------#

# Simulate the population dynamics of the 11 empirical multilayer networks
# under increased hurricane intensity, the second scenario of the Methods.
# Every island keeps its empirical hurricane frequency lambda, and a shift s
# moves a share of the probability mass of its PMF from categories 1 to 3 to
# categories 4 and 5. The shift sweeps from 0 to 1 in steps of 0.05, so each
# island runs 21 scenarios of 1,000 replicates of 300 years, 231 tasks in
# all. The extinction curves against the shift are Fig. 3d (survival_curves.R
# and figure_3.R) and give the species fragility to hurricane intensity of
# the pSEM in Fig. 4 (calculate_auc_extinction.R). The shift of 0 reproduces
# the empirical regime of each island, the baseline of the null-model
# comparison in Extended Data Fig. 1 (null_model_comparison.R).
#
# Inputs:  src/multi_pol/B_<island>.csv and src/multi_sd/S_<island>.csv, the
#          incidence matrices of the two layers, built by
#          scripts/01_network_preparation/build_networks_multilayer.R
#          src/empirical_lambdas.csv and src/empirical_pmfs.csv, written by
#          empirical_lambda_pmf.R
# Outputs: data/main_sims_pmf/<island>_multilayer_lambda<lambda>_pmf<s>_tyears300.csv,
#          one file per island and shift s, read by survival_curves.R and,
#          for s = 0, by null_model_comparison.R
# Helpers: src/main_functions.jl (population dynamics, hurricanes and the
#          simulation wrapper) and src/aux_functions.jl (networks, regimes
#          and the interaction matrix). main_simulations_empirical.jl runs
#          the companion sweep over hurricane frequency.

using Distributed

# Worker processes; pmap gives each one scenario at a time
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
    #-2. Empirical regime of each island---------------------------#
    #--------------------------------------------------------------#

    # Frequency and PMF of each island, keyed by the island code
    regimes = empirical_regimes()

    #--------------------------------------------------------------#
    #-3. Shift of probability mass toward categories 4 and 5-------#
    #--------------------------------------------------------------#

    # Shift probability mass from categories 1 to 3 to categories 4 and 5, the
    # increased hurricane intensity of the Methods. Within each group the
    # categories keep their empirical proportions, so the shift raises the
    # intensity of hurricanes without reordering the categories. The shift
    # cannot exceed the mass of categories 1 to 3, so a shift of 1 moves all
    # of it.
    # Arguments: pmf maps each category of one island to its probability and
    # damage fraction, as in empirical_regimes; shift is the share of
    # probability mass moved. Returns a Dict of the same form with the shifted
    # probabilities.
    function shift_pmf(; pmf::Dict, shift::Float64)
        low = ("Cat_1", "Cat_2", "Cat_3")
        p_low = sum(pmf[cat].prob for cat in low if haskey(pmf, cat); init=0.0)
        p_high = sum(pmf[cat].prob for cat in ("Cat_4", "Cat_5") if haskey(pmf, cat); init=0.0)
        moved = min(shift, p_low)

        # Categories 1 to 3 lose the shifted mass, and categories 4 and 5
        # receive it, each in proportion to their probabilities
        scale_low = p_low > 0.0 ? (p_low - moved) / p_low : 0.0
        scale_high = p_high > 0.0 ? (p_high + moved) / p_high : 0.0

        return Dict(cat => (prob=v.prob * (cat in low ? scale_low : scale_high), damage=v.damage)
                    for (cat, v) in pmf)
    end

    #--------------------------------------------------------------#
    #-4. One scenario per island and shift-------------------------#
    #--------------------------------------------------------------#

    # Intensity shifts of the sweep. Each step moves another 5% of the
    # probability mass to categories 4 and 5, and the shift of 0 reproduces the
    # empirical regime of each island.
    shifts = collect(0.0:0.05:1.0)

    # Each island meets its own empirical frequency and its shifted PMF
    c_list = vec([Dict(
        :net => net,
        :lambda => regimes[net.network].lambda,
        :intensity_scenario => shift_pmf(pmf=regimes[net.network].pmf, shift=shift),
        :shift => shift,
        :m => 0.5,
        :d => 0.5,
        :c => -0.1,
        :ht => 0.5,
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
