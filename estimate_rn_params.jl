using Catalyst
using DifferentialEquations
using Zygote
using SciMLSensitivity
using Optimization
using OptimizationOptimisers
using ForwardDiff
using Plots

include("construct_rn.jl")

function sample_rn(rn::ReactionSystem,
                   topology::String,
                   k_ons::Vector,
                   t_span::Tuple;
                   num_samples::Int64=100,
                   noise::Float64=0.1)::Tuple{Vector, Array, ODEProblem}
    # Start with concentration 1 for all monomers, whose names have 5 characters,
    # e.g. "X1(t)"
    u0 = [length(string(species)) == 5 for species in species(rn)]
    # The time span should start from a small, but non-zero time in order to avoid 
    # taking the log of zero.
    # Obtain a log10-spaced vector of times at which to sample the ODE solution.
    # The log10-scale enables us to observe trapping at early simulation times.
    sample_times = 10 .^ range(log10(t_span[1]), log10(t_span[2]), length=num_samples)

    params = get_rate_constants_from_k_ons(k_ons, topology)
    prob = ODEProblem(rn, u0, t_span, params)
    sol_real = solve(prob, Tsit5(); tstops = sample_times)
    sample_vals = Array(sol_real(sample_times))

    # Add noise to the ODE solution by uniformly sampling a coefficient for the
    # observed value at each time step from the interval [1-`noise`/2, 1+`noise`/2].
    sample_vals .*= (1 .+ noise * rand(Float64, size(sample_vals)) .- (noise / 2))

    return sample_times, sample_vals, prob
end

function detect_traps(t::Vector, u::Vector)
    # Approximate the time derivative as a simple finite difference
    finite_difference = diff(u) ./ diff(t)
    trap_bounds = []
    in_trap = false
    trap_start = 0.0
    trap_end = 0.0
    threshold = std(u)
    for i in eachindex(finite_difference[1:end])
        if finite_difference[i] < threshold && !in_trap
            in_trap = true
            trap_start = t[i]
        end
        if finite_difference[i] > threshold && in_trap
            in_trap = false
            trap_end = t[i]
            push!(trap_bounds, (trap_start, trap_end))
        end
    end
    if in_trap
        push!(trap_bounds, (trap_start, t[end]))
    end

    print(trap_bounds)
end

function optimise_p(p_init::Vector, 
                    topology::String,
                    t_span::Tuple,
                    sample_times::Vector{Float64}, 
                    sample_vals::Array{Float64, 2}, 
                    prob::ODEProblem)::Vector{Float64}
    """
    `topology` should be "fc" for fully connected topology or "rg" for ring topology
    """
    function loss(p, _)
        newtimes = filter(<=(t_span[2]), sample_times)
        params = get_rate_constants_from_k_ons(p, topology)
        newprob = remake(prob; tspan=t_span, p=params)
        sol = Array(solve(newprob, TRBDF2(); saveat = newtimes))

        # Optimize based only on end product
        #loss = sum(abs2, sol[7, :] .- sample_vals[7, 1:size(sol,2)]) 
        # Optimize based on all species --- this unevenly weights moderate levels of assembly with weight (n choose k)
        loss = sum(abs2, sol .- sample_vals[:, 1:size(sol,2)])

        return loss, sol
    end

    # optimize for the parameters that minimize the loss
    optf = OptimizationFunction(loss, Optimization.AutoZygote())
    optprob = OptimizationProblem(optf, p_init)
    sol = solve(optprob, OptimizationOptimisers.ADAM(10); maxiters = 200)

    # Return the parameter estimates
    return sol.u
end

function plot_rn(rn::ReactionSystem,
    topology::String,
    p::Vector;
    plot_type::String="line",
    t_span::Tuple=(1e-4, 1e5),
    noise::Float64=0.,
    num_samples::Int64=100,
    linesize::Int64=4)
    sample_times, sample_vals, _ = 
    sample_rn(rn, topology, p, t_span, noise=noise, num_samples=num_samples)

    if plot_type == "line"
        plot(log10.(sample_times), sample_vals[1, :], lw=linesize, label="Monomer")
        plot!(log10.(sample_times), sample_vals[3, :], lw=linesize, label="Dimer")
        plot!(log10.(sample_times), sample_vals[7, :], lw=linesize, label="Trimer")
    elseif plot_type == "dot"
        scatter(log10.(sample_times), sample_vals[1, :], lw=linesize, label="Monomer")
        scatter!(log10.(sample_times), sample_vals[3, :], lw=linesize, label="Dimer")
        scatter!(log10.(sample_times), sample_vals[7, :], lw=linesize, label="Trimer")
    end
    plot!(xlabel="log10(t)", ylabel="Relative yield")
end