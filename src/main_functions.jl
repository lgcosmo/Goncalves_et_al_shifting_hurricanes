#--------------------------------------------------------------#
#-Population dynamics, hurricanes and the simulation wrapper---#
#--------------------------------------------------------------#

# Functions that simulate the population dynamics of one multilayer
# bird-pollination and seed-dispersal network under a hurricane regime, as
# in the Methods. Species abundances follow ordinary differential equations
# with an intrinsic growth rate, intra-guild competition, and mutualistic
# benefits from pollination and seed dispersal. The benefits saturate
# through a Type II functional response and are diluted by the number of
# partners through a generalist trade-off. Every network first grows to its
# steady state. Then hurricanes arrive as a Poisson process with frequency
# lambda, each removes a damage fraction drawn from the island's PMF of
# categories from every species, and a species whose abundance falls below
# its extinction threshold goes extinct. The replicates of one scenario form
# one simulation file, the raw material of the extinction curves (Fig. 3c
# and 3d), of species fragility (Fig. 4) and of the null-model comparison
# (Extended Data Fig. 1).
#
# Inputs:  one network from load_network and the hurricane regime of one
#          scenario, passed by the simulation scripts
# Outputs: data/<out_dir>/<network>_multilayer_lambda<lambda>_pmf<i>_tyears<t_years>.csv,
#          one row per species and replicate, read by survival_curves.R,
#          process_null_simulations.R and null_model_comparison.R
# Used by: main_simulations_empirical.jl (frequency sweep),
#          main_simulations_empirical_pmf.jl (intensity sweep) and
#          main_simulations_null_decoupling.jl (null networks) in
#          scripts/03_simulations, all through net_multisim_typeII
# Helpers: load_network and yij_multi! in src/aux_functions.jl. The wrapper
#          also needs CSV and DrWatson, which the scripts load.

using DifferentialEquations
using Distributions
using DataFrames
using Random
using LinearAlgebra
using Parameters

#--------------------------------------------------------------#
#-Population dynamics------------------------------------------#
#--------------------------------------------------------------#

# The population dynamics model of the Methods, in the form the ODE solver
# calls. Species i grows at its intrinsic growth rate r_i, is limited
# by its own abundance, loses abundance to its intra-guild competitors, and
# gains mutualistic benefits from its pollination and seed-dispersal
# partners. The benefits saturate through a Type II functional response with
# handling time h, so that abundant partners bring diminishing returns.
# Arguments: dN receives the derivatives; N holds the current abundances; p
# is the tuple (Y_mut, Y_comp, r, h, mut_raw, comp_sum, M) with the
# mutualistic and the competitive part of the interaction matrix, the growth
# rates, the handling times and three work vectors that spare allocations
# during the integration; t is the time, which the solver passes and the
# dynamics ignore. Fills dN in place and returns it.
function eco_dynamics_typeII!(dN::Vector{Float64}, N::Vector{Float64}, p::Tuple, t::Float64)
    Y_mut::Matrix{Float64}, Y_comp::Matrix{Float64}, r::Vector{Float64}, h::Vector{Float64}, mut_raw::Vector{Float64}, comp_sum::Vector{Float64}, M::Vector{Float64} = p

    # Total benefit received from partners and total effect of competitors
    mul!(mut_raw, Y_mut, N)
    mul!(comp_sum, Y_comp, N)

    # The benefit saturates with handling time h (Type II functional response)
    @inbounds @. M = mut_raw / (1.0 + h * mut_raw)

    # Growth of each species, limited by its own abundance and its competitors
    # and raised by its saturated benefit
    @inbounds for i in eachindex(N)
        dN[i] = N[i] * (r[i] - N[i] + comp_sum[i] + M[i])
    end

    return dN
end

#--------------------------------------------------------------#
#-Steady state, hurricane damage, extinction and recovery------#
#--------------------------------------------------------------#

# Condition of the callback that ends the growth to the steady state. The
# network has reached its steady state when the norm of the abundance
# changes falls below 1e-3. Arguments: N holds the current abundances; t is
# the time; integrator is the solver object, whose field p holds the
# parameters of the dynamics. Returns the norm minus 1e-3, which crosses
# zero at the steady state.
function steady_state_condition_typeII(N, t, integrator)
    dN = similar(N)
    eco_dynamics_typeII!(dN, N, integrator.p, t)
    return norm(dN) - 1e-3
end

# Action of the callback. It ends the integration at the steady state, so
# that the hurricane regime starts from the steady-state abundances.
# Arguments: integrator is the solver object.
function stop_integration!(integrator)
    DifferentialEquations.terminate!(integrator)
end

# Sampler of the damage fraction of one hurricane. The category of each
# hurricane follows the island's probability mass function (PMF) of
# Saffir-Simpson categories, the hurricane intensity of the Methods, and each
# category carries its own damage fraction. The sampler is built once per
# replicate, and each hurricane then draws one category from it.
# Arguments: pmf maps each category name ("Cat_1" to "Cat_5") to a NamedTuple
# with the probability of the category (prob) and its damage fraction
# (damage). Returns a NamedTuple with the categorical distribution of the
# categories (categories) and their damage fractions in the same order
# (damages), so that damages[rand(categories)] is the damage fraction of one
# hurricane.
function damage_sampler(; pmf::Dict)
    names = collect(keys(pmf))
    return (categories=Categorical([pmf[cat].prob for cat in names]),
            damages=[pmf[cat].damage for cat in names])
end

# Time at which each species went extinct, counted from the start of the
# hurricane regime. A hurricane that pushes a species below its extinction
# threshold sets its abundance to zero, so the first recorded abundance
# below threshold marks the extinction. survival_curves.R counts the
# extinctions of the frequency and intensity sweeps from these times.
# Arguments: R holds the recorded abundances (species in rows, times in
# columns); all_t holds the recorded times; start_time is the time at which
# the network reached its steady state and the hurricane regime began;
# threshold is the abundance that counts as zero (1e-9 by default). Returns
# one time per species, and Inf for the species that survived.
function get_extinction_times(; R::Matrix{Float64}, all_t::Vector{Float64}, start_time::Float64, threshold::Float64=1e-9)
    nsp = size(R, 1)
    relative_extinction_times = fill(Inf, nsp)

    @inbounds for i in 1:nsp
        ext_idx = findfirst(R[i, :] .< threshold)
        if !isnothing(ext_idx)
            absolute_ext_time = all_t[ext_idx]
            relative_time = absolute_ext_time - start_time
            relative_extinction_times[i] = max(0.0, relative_time)
        end
    end
    return relative_extinction_times
end

# Mean recovery rate of each species over the hurricanes of a replicate. The
# recovery rate after one hurricane is the change in log abundance over the
# recovery window, divided by the length of the window. The window starts
# from the abundance right after the hurricane, which the record holds as
# the last point at the time of the hurricane.
# Arguments: R and all_t are as in get_extinction_times; perturbation_times
# holds the times of the hurricanes; recovery_window is the length of the
# window in years (5 by default). Returns one rate per species, averaged over
# the hurricanes it survived with a complete window, and NaN for a species
# with no such hurricane. The simulation files carry this rate, but no later
# script reads it.
function get_avg_recovery_rates(; R::Matrix{Float64}, all_t::Vector{Float64}, perturbation_times::Vector{Float64}, recovery_window::Float64=5.0)
    nsp = size(R, 1)
    avg_recovery_rates = fill(NaN, nsp)

    for i in 1:nsp
        species_rates = Float64[]
        for event_time in perturbation_times
            # The record holds the time of a hurricane twice, before and after
            # its damage, and the window starts after the damage
            eq_inds = findall(x -> x == event_time, all_t)
            if !isempty(eq_inds)
                start_idx = eq_inds[end]
            else
                start_idx = findfirst(all_t .> event_time)
            end

            end_idx = findfirst(all_t .>= (event_time + recovery_window))

            if !isnothing(start_idx) && !isnothing(end_idx)
                N_start = R[i, start_idx]
                N_end = R[i, end_idx]

                # Only a species alive at both ends of the window has a rate
                if N_start > 0.0 && N_end > 0.0
                    rate = (log(N_end) - log(N_start)) / recovery_window
                    push!(species_rates, rate)
                end
            end
        end

        if !isempty(species_rates)
            avg_recovery_rates[i] = mean(species_rates)
        end
    end
    return avg_recovery_rates
end

#--------------------------------------------------------------#
#-One replicate of one network---------------------------------#
#--------------------------------------------------------------#

# Simulate one replicate of one network under one hurricane regime, the
# event-driven simulation of the Methods. The replicate has two phases.
# First, the network grows from random initial abundances to its steady
# state, without hurricanes. Second, hurricanes strike at random times for
# t_years years. The waiting time between hurricanes is exponential with
# mean 1/lambda, as in a Poisson process, and each hurricane removes a
# damage fraction drawn from the intensity PMF from every species. A species
# whose abundance falls below its extinction threshold, 1% of its
# steady-state abundance, goes extinct.
# Arguments: net is the network from load_network, with its species, their
# category and the binary matrices of the two layers; lambda is the
# hurricane frequency; intensity_pmf maps each category to its probability
# and damage fraction, as in damage_sampler; intensity labels the scenario in
# the output (the shift of an intensity scenario, or the expected damage
# fraction of the PMF); t_years is the length of the hurricane regime; m, d
# and c are the base strengths of pollination, seed dispersal and intra-guild
# competition; ht is the handling time of the Type II functional response;
# sim is the number of the replicate, written to the output and not used as
# a seed.
# Returns a DataFrame with one row per species. It holds the intrinsic growth
# rate, the abundance at the end of the record, whether and when the species
# went extinct, the time at which the network reached its steady state, the
# mean recovery rate, the parameters, the hurricane frequency and the
# intensity label.
function net_ecoperturb_typeII(; net::NamedTuple, lambda::Float64, intensity_pmf::Dict, intensity::Float64, t_years::Int64, m::Float64, c::Float64, d::Float64, ht::Float64, sim::Int64)

    nt = length(net.species)
    Y = zeros(Float64, nt, nt)
    # Intrinsic growth rates between 0 and 0.5 and initial abundances between
    # 0.1 and 1, drawn at random for each replicate
    r = rand(0.0:0.01:0.5, nt)
    N0 = rand(0.1:0.01:1.0, nt)
    h = fill(ht, nt)

    # Work vectors of the dynamics, allocated once per replicate
    mut_raw = zeros(Float64, nt)
    comp_sum = zeros(Float64, nt)
    M_vec = zeros(Float64, nt)

    # Interaction matrix of the multilayer network, split into its mutualistic
    # (positive) and competitive (negative) parts
    yij_multi!(B=net.B, S=net.S, Y=Y, m=m, c=c, d=d)

    Y_mut = Y .* Float64.(Y .> 0.0)
    Y_comp = Y .* Float64.(Y .<= 0.0)
    p_dynamics = (Y_mut, Y_comp, r, h, mut_raw, comp_sum, M_vec)

    # First phase. The network grows to its steady state within at most 1,000
    # years, without hurricanes.
    tspan_initial = (0.0, 1000.0)
    model_template = ODEProblem(eco_dynamics_typeII!, N0, tspan_initial, p_dynamics)
    cb = ContinuousCallback(steady_state_condition_typeII, stop_integration!)
    sol_initial = solve(model_template, callback=cb)

    # Second phase. Hurricanes strike for t_years years, starting from the
    # steady state.
    current_time = sol_initial.t[end]
    t_end_final = current_time + t_years
    current_state = sol_initial.u[end]
    neq_initial = sol_initial.u[end]
    epsilon = 0.01 # Extinction threshold of each species, 1% of its steady-state abundance
    extinction_thresholds = neq_initial .* epsilon
    runs = [sol_initial]
    perturbation_times = Float64[]
    event_dist = Exponential(1.0 / lambda)
    sampler = damage_sampler(pmf=intensity_pmf)

    # Hurricanes arrive one at a time, with exponential waiting times. The
    # regime ends when the next hurricane would fall after t_years, so the
    # record of the replicate ends at the last hurricane, before its damage.
    while current_time < t_end_final
        time_to_next_hurricane = rand(event_dist)
        event_time = current_time + time_to_next_hurricane
        if event_time > t_end_final
            break
        end

        # Abundances recover until the hurricane strikes
        tspan_wait = (current_time, event_time)
        model_wait = remake(model_template; u0=current_state, tspan=tspan_wait)
        sol_wait = solve(model_wait)
        push!(runs, sol_wait)

        # The hurricane removes the damage fraction of its category from every
        # species, and the species that fall below their extinction threshold
        # go extinct
        state_before_hit = sol_wait.u[end]
        damage_fraction = sampler.damages[rand(sampler.categories)]

        state_after_hit = state_before_hit .* (1.0 - damage_fraction)
        state_after_hit[state_after_hit .< extinction_thresholds] .= 0.0

        current_state = state_after_hit
        current_time = event_time
        push!(perturbation_times, current_time)
    end

    # Join the records of both phases and measure extinctions and recovery
    perturbation_start_time = sol_initial.t[end]
    all_t = copy(runs[1].t)
    all_u = copy(runs[1].u)
    for i in 2:length(runs)
        append!(all_t, runs[i].t)
        append!(all_u, runs[i].u)
    end

    val_matrix = hcat(all_u...)
    ext_times = get_extinction_times(R=val_matrix, all_t=all_t, start_time=perturbation_start_time)
    avg_recovery = get_avg_recovery_rates(R=val_matrix, all_t=all_t, perturbation_times=perturbation_times)

    df = DataFrame(
        network = net.network,
        sim_id = sim,
        species_id = net.species,
        type = net.category,
        r = r,
        final_abundance = round.(val_matrix[:, end], digits=4),
        is_extinct = isfinite.(ext_times),
        time_to_extinction = round.(ext_times, digits=4),
        initial_teq = round.(perturbation_start_time, digits=4),
        avg_recovery_rate = round.(avg_recovery, digits=4),
        m = m, d = d, c = c,
        intensity = intensity,
        lambda = lambda
    )

    return df
end

#--------------------------------------------------------------#
#-All replicates of one scenario-------------------------------#
#--------------------------------------------------------------#

# Run every replicate of one scenario and write them to one simulation file,
# which survival_curves.R and process_null_simulations.R read. The
# simulation scripts call this function through pmap, one scenario per task,
# and the replicates of a scenario run one after the other on its worker.
# Arguments: p is the dictionary of one scenario, built by dict_list or by
# hand. It holds nsim, the number of replicates; net, the network from
# load_network; m, c, d and ht, the parameters of the dynamics; lambda, the
# hurricane frequency; intensity_scenario, the PMF of categories; and
# t_years, the length of the hurricane regime. It may also hold shift, the
# intensity shift of the scenario, and out_dir, the folder under data/ that
# receives the file ("sims" by default).
# Writes data/<out_dir>/<network>_multilayer_lambda<lambda>_pmf<i>_tyears<t_years>.csv,
# where i is the intensity label, the shift of an intensity scenario or,
# without a shift, the expected damage fraction of the PMF. Returns nothing.
function net_multisim_typeII(p::Dict)
    @unpack nsim, net, m, c, d, ht, lambda, intensity_scenario, t_years = p

    shift = haskey(p, :shift) ? p[:shift] : NaN
    out_dir = haskey(p, :out_dir) ? p[:out_dir] : "sims"

    # Intensity label of the scenario, written to the output and the file name
    intensity = isnan(shift) ? round(sum(v.prob * v.damage for v in values(intensity_scenario)), digits=2) : shift

    results_list = Vector{DataFrame}(undef, nsim)
    println("Starting $(nsim) simulations for network: $(net.network) with PMF: $(intensity)...")

    for n in 1:nsim
        results_list[n] = net_ecoperturb_typeII(
            net = net, lambda = lambda, intensity_pmf = intensity_scenario,
            intensity = intensity, t_years = t_years, m = m, c = c, d = d, ht = ht, sim = n
        )
    end
    println("...simulations complete.")

    final_results = vcat(results_list...)
    filename = "$(net.network)_multilayer_lambda$(lambda)_pmf$(intensity)_tyears$(t_years).csv"
    mkpath(datadir(out_dir))
    CSV.write(datadir(out_dir, filename), final_results)
    println("Results saved to $(joinpath(out_dir, filename))")
end
