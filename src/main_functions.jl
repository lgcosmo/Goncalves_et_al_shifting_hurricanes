#--------------------------------------------------------------#
#-Population dynamics, hurricanes and the simulation wrapper---#
#--------------------------------------------------------------#

# Population dynamics of one multilayer network (pollination and seed
# dispersal) under a hurricane regime. Abundances follow ODEs with intrinsic
# growth, intra-guild competition and mutualistic benefits that saturate
# through a Type II functional response and are diluted by the number of
# partners (generalist trade-off). Each network first grows to its steady
# state. Then hurricanes arrive as a Poisson process with frequency lambda,
# each removes a damage fraction drawn from the PMF of categories of the
# island, and species below their extinction threshold go extinct.
#
# Inputs:  one network from load_network (src/aux_functions.jl) and the
#          hurricane regime of one scenario
# Outputs: one file per scenario in data/<out_dir>/, named
#          <network>_multilayer_lambda<lambda>_pmf<i>_tyears<t_years>.csv,
#          with one row per species and replicate

using DifferentialEquations
using Distributions
using DataFrames
using Random
using LinearAlgebra
using Parameters

#--------------------------------------------------------------#
#-Population dynamics------------------------------------------#
#--------------------------------------------------------------#

# Population dynamics in the form called by the ODE solver. Species i grows
# at rate r_i, is limited by its own abundance and by its intra-guild
# competitors, and gains benefits from its pollination and seed-dispersal
# partners. p is the tuple (Y_mut, Y_comp, r, h, mut_raw, comp_sum, M):
# mutualistic and competitive parts of the interaction matrix, growth rates,
# handling times and three preallocated work vectors.
function eco_dynamics_typeII!(dN::Vector{Float64}, N::Vector{Float64}, p::Tuple, t::Float64)
    Y_mut::Matrix{Float64}, Y_comp::Matrix{Float64}, r::Vector{Float64}, h::Vector{Float64}, mut_raw::Vector{Float64}, comp_sum::Vector{Float64}, M::Vector{Float64} = p

    # Total benefit received from partners and total effect of competitors
    mul!(mut_raw, Y_mut, N)
    mul!(comp_sum, Y_comp, N)

    # The benefit saturates with handling time h (Type II functional response)
    @inbounds @. M = mut_raw / (1.0 + h * mut_raw)

    # Growth of each species: r - N + competition + saturated benefit
    @inbounds for i in eachindex(N)
        dN[i] = N[i] * (r[i] - N[i] + comp_sum[i] + M[i])
    end

    return dN
end

#--------------------------------------------------------------#
#-Steady state, hurricane damage, extinction and recovery------#
#--------------------------------------------------------------#

# Callback condition that ends the growth to the steady state: the norm of
# the abundance changes falls below 1e-3.
function steady_state_condition_typeII(N, t, integrator)
    dN = similar(N)
    eco_dynamics_typeII!(dN, N, integrator.p, t)
    return norm(dN) - 1e-3
end

# Callback action: stop the integration at the steady state.
function stop_integration!(integrator)
    DifferentialEquations.terminate!(integrator)
end

# Sampler of the damage fraction of one hurricane. pmf maps each category
# ("Cat_1" to "Cat_5") to its probability (prob) and damage fraction
# (damage). Returns the categorical distribution of the categories and their
# damage fractions in the same order, so damages[rand(categories)] is the
# damage of one hurricane.
function damage_sampler(; pmf::Dict)
    names = collect(keys(pmf))
    return (categories=Categorical([pmf[cat].prob for cat in names]),
            damages=[pmf[cat].damage for cat in names])
end

# Time at which each species went extinct, counted from the start of the
# hurricane regime (start_time). A hurricane that pushes a species below its
# extinction threshold sets its abundance to zero, so the first abundance
# below threshold (1e-9) is the extinction. R has species in rows and
# recorded times in columns. Returns Inf for the species that survived.
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
# recovery window (5 years), divided by the length of the window. Returns
# NaN for a species with no complete window.
function get_avg_recovery_rates(; R::Matrix{Float64}, all_t::Vector{Float64}, perturbation_times::Vector{Float64}, recovery_window::Float64=5.0)
    nsp = size(R, 1)
    avg_recovery_rates = fill(NaN, nsp)

    for i in 1:nsp
        species_rates = Float64[]
        for event_time in perturbation_times
            # The time of a hurricane is recorded twice, before and after the
            # damage. The window starts after the damage.
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

# Simulate one replicate of one network under one hurricane regime. First
# the network grows from random initial abundances to its steady state,
# without hurricanes. Then hurricanes arrive for t_years years, with
# exponential waiting times of mean 1/lambda (Poisson process). Each
# hurricane removes a damage fraction drawn from intensity_pmf from every
# species, and a species below its extinction threshold (1% of its
# steady-state abundance) goes extinct.
# m, d and c are the base strengths of pollination, seed dispersal and
# intra-guild competition, ht the handling time, delta the exponent of the
# generalist trade-off and tradeoff the degree it uses ("layer" or
# "multilayer", see yij_multi!). intensity and sim are only written to the
# output. Returns a DataFrame with one row per species.
function net_ecoperturb_typeII(; net::NamedTuple, lambda::Float64, intensity_pmf::Dict, intensity::Float64, t_years::Int64, m::Float64, c::Float64, d::Float64, ht::Float64, delta::Float64, tradeoff::String, sim::Int64)

    nt = length(net.species)
    Y = zeros(Float64, nt, nt)
    # Random growth rates (0 to 0.5) and initial abundances (0.1 to 1)
    r = rand(0.0:0.01:0.5, nt)
    N0 = rand(0.1:0.01:1.0, nt)
    h = fill(ht, nt)

    # Preallocated work vectors of the dynamics
    mut_raw = zeros(Float64, nt)
    comp_sum = zeros(Float64, nt)
    M_vec = zeros(Float64, nt)

    # Interaction matrix, split into its mutualistic (positive) and
    # competitive (negative) parts
    yij_multi!(B=net.B, S=net.S, Y=Y, m=m, c=c, d=d, delta=delta, tradeoff=tradeoff)

    Y_mut = Y .* Float64.(Y .> 0.0)
    Y_comp = Y .* Float64.(Y .<= 0.0)
    p_dynamics = (Y_mut, Y_comp, r, h, mut_raw, comp_sum, M_vec)

    # First phase: growth to the steady state, without hurricanes (at most
    # 1,000 years)
    tspan_initial = (0.0, 1000.0)
    model_template = ODEProblem(eco_dynamics_typeII!, N0, tspan_initial, p_dynamics)
    cb = ContinuousCallback(steady_state_condition_typeII, stop_integration!)
    sol_initial = solve(model_template, callback=cb)

    # Second phase: hurricanes for t_years years, from the steady state
    current_time = sol_initial.t[end]
    t_end_final = current_time + t_years
    current_state = sol_initial.u[end]
    neq_initial = sol_initial.u[end]
    epsilon = 0.01 # Extinction threshold: 1% of the steady-state abundance
    extinction_thresholds = neq_initial .* epsilon
    runs = [sol_initial]
    perturbation_times = Float64[]
    event_dist = Exponential(1.0 / lambda)
    sampler = damage_sampler(pmf=intensity_pmf)

    # Hurricanes arrive one at a time with exponential waiting times, until
    # the next one would fall after t_years. The record ends at the last
    # hurricane, before its damage.
    while current_time < t_end_final
        time_to_next_hurricane = rand(event_dist)
        event_time = current_time + time_to_next_hurricane
        if event_time > t_end_final
            break
        end

        # Integrate the dynamics until the next hurricane
        tspan_wait = (current_time, event_time)
        model_wait = remake(model_template; u0=current_state, tspan=tspan_wait)
        sol_wait = solve(model_wait)
        push!(runs, sol_wait)

        # The hurricane removes its damage fraction from every species.
        # Species below their extinction threshold go extinct.
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
        delta = delta, tradeoff = tradeoff,
        intensity = intensity,
        lambda = lambda
    )

    return df
end

#--------------------------------------------------------------#
#-All replicates of one scenario-------------------------------#
#--------------------------------------------------------------#

# Run all replicates of one scenario and save them to one simulation file.
# p is the dictionary of the scenario, with nsim (number of replicates), net,
# the parameters of the dynamics (m, c, d, ht, delta, tradeoff), lambda,
# intensity_scenario (PMF of categories) and t_years. Optionally shift (the
# intensity shift) and out_dir (folder under data/, "sims" by default).
# The file name has only the network and the hurricane regime, so scenarios
# that differ in the parameters of the dynamics need their own out_dir.
function net_multisim_typeII(p::Dict)
    @unpack nsim, net, m, c, d, ht, delta, tradeoff, lambda, intensity_scenario, t_years = p

    shift = haskey(p, :shift) ? p[:shift] : NaN
    out_dir = haskey(p, :out_dir) ? p[:out_dir] : "sims"

    # Intensity label: the shift or, without a shift, the expected damage
    # fraction of the PMF
    intensity = isnan(shift) ? round(sum(v.prob * v.damage for v in values(intensity_scenario)), digits=2) : shift

    results_list = Vector{DataFrame}(undef, nsim)
    println("Starting $(nsim) simulations for network: $(net.network) with PMF: $(intensity)...")

    for n in 1:nsim
        results_list[n] = net_ecoperturb_typeII(
            net = net, lambda = lambda, intensity_pmf = intensity_scenario,
            intensity = intensity, t_years = t_years, m = m, c = c, d = d, ht = ht,
            delta = delta, tradeoff = tradeoff, sim = n
        )
    end
    println("...simulations complete.")

    final_results = vcat(results_list...)
    filename = "$(net.network)_multilayer_lambda$(lambda)_pmf$(intensity)_tyears$(t_years).csv"
    mkpath(datadir(out_dir))
    CSV.write(datadir(out_dir, filename), final_results)
    println("Results saved to $(joinpath(out_dir, filename))")
end
