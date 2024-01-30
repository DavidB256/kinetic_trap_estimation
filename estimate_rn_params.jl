using Catalyst
using DifferentialEquations
using Optimization
using OptimizationOptimisers
using ForwardDiff
using PreallocationTools

# These imports may be unncessary
using Zygote
using SciMLSensitivity
using OptimizationOptimJL

include("construct_rn.jl")

function log10_range(start, ending, length::Int64=100)::Vector
    return 10 .^ range(log10(start), log10(ending), length=length)
end

function sample_rn(rn::ReactionSystem,
                   topology::String,
                   k_ons::Vector,
                   sample_times::Vector;
                   noise::Float64=0.,
                   printer::Bool=false)

    if printer
        print("True\t")
        for k_on in k_ons
            print(k_on)
            print("\t")
        end
        println()
    end

    # Start with concentration 1 for all monomers, whose names have 5 characters,
    # e.g. "X1(t)"
    u0 = [length(string(species)) == 5 for species in species(rn)]
    # The time span should start from a small, but non-zero time in order to avoid 
    # taking the log of zero.
    # Obtain a log10-spaced vector of times at which to sample the ODE solution.
    # The log10-scale enables us to observe trapping at early simulation times.
    t_span = (sample_times[1], sample_times[end])

    params = get_rate_constants_from_k_ons(k_ons, topology)
    prob = ODEProblem(rn, u0, t_span, params)
    sol_real = solve(prob, TRBDF2(); tstops = sample_times)
    sample_vals = Array(sol_real(sample_times))

    # Add noise to the ODE solution by uniformly sampling a coefficient for the
    # observed value at each time step from the interval [1-`noise`/2, 1+`noise`/2].
    sample_vals .*= (1 .+ noise .* rand(Float64, size(sample_vals)) .- (noise / 2))

    return sample_vals, prob
end

# Custom type used to record the height and horizontal bounds of kinetic traps
mutable struct kinetic_trap
    y::Float64
    t_start::Float64
    t_end::Float64
end

"""
    detect_traps(t, u)

Given equal-length vectors of times `t` concentrations `u`, return a vector
of `trap` objects that record the height and horizontal bounds of kinetic traps,
i.e. regions in which the finite difference of `u` w.r.t. `log10.(t)` is less than
`threshold=0.05`.
"""
function detect_traps(t::Vector, u::Vector)::Vector{kinetic_trap}
    # Approximate the time derivative as a simple finite difference
    
    finite_difference = diff(u) ./ diff(log10.(t))
    kinetic_traps::Vector{kinetic_trap} = []
    trap_start_index::Float64 = 0.
    in_trap = false # prevent counting beginning of a reaction as a trap
    has_exceeded_threshold = false
    threshold = 0.1

    current_trap = kinetic_trap(0., 0., 0.)
    for i in eachindex(finite_difference[1:end])
        if (finite_difference[i] < threshold) && !in_trap && has_exceeded_threshold
            in_trap = true
            trap_start_index = i
            current_trap.t_start = t[i]
        end
        if (finite_difference[i] > threshold) && in_trap
            in_trap = false
            current_trap.t_end = t[i]
            current_trap.y = u[Int(round((i + trap_start_index) / 2))]
            push!(kinetic_traps, current_trap)
            current_trap = kinetic_trap(0., 0., 0.) # Clean out fields, just to be safe
        end
        if finite_difference[i] > threshold
            has_exceeded_threshold = true
        end
    end

    return kinetic_traps
end

function optimise_p(p_init::Vector, 
                    topology::String,
                    sample_times::Vector{Float64}, 
                    sample_vals::Array{Float64, 2}, 
                    prob::ODEProblem,
                    autodiff_routine::String="reverse",
                    experiment_type::String="end_product_yield")::Vector{Float64}
    # `topology` should be "fc" for fully connected topology or "rg" for ring topology

    if experiment_type == "kinetic_traps"
        # Ensure 100 sample times for trap detection
        sample_times = log10_range(sample_times[1], sample_times[end])
        # sample_vals should be noise-free
        kinetic_traps = detect_traps(sample_times, sample_vals[end, :])
    end

    function loss(p::AbstractVector{T}, _) where T
        params = get_rate_constants_from_k_ons(p, topology)

        params = T[params[1],
                   params[2],
                   params[3],
                   params[4],
                   params[5],
                   params[6],
                   params[7],
                   params[8]]

        #println("params in loss func: $(params)")

        newprob = remake(prob; 
                         p=params, 
                         tspan=T[sample_times[1], sample_times[end]])
        # TRBDF2 is a stiff solver that is warranted by the log-time-scale of kinetic
        # kinetic traps, there are situations in which Tsit5 can be used instead
        sol = Array(solve(newprob, TRBDF2(); saveat=sample_times, maxiters=1e7))

        # sol_vals = zeros(T, size(sol, 2))
        # sample_vals = zeros(T, size(sol, 2))

        # for i = 1:size(sol, 2)
        #     sol_vals[i] = sol[7, i]
        #     sample_vals[i] = sample_vals[7, i]
        # end

        if experiment_type == "end_product_yield"
            loss = sum((sol[end, :] .- sample_vals[end, :]) .^ 2)
        elseif experiment_type == "kinetic_traps"
            loss = 0
            for kt in kinetic_traps
                idxs = findall(sample_times .>= kt.t_start .&& sample_times .<= kt.t_end)

                loss += sum((kt.y .- sol[end, idxs]) .^ 2)
            end
        else
            error("Invalid experiment type: $(experiment_type)")
        end
        
        return loss, sol
    end

    # optimize for the parameters that minimize the loss
    if autodiff_routine in ["rev", "reverse"]
        optf = OptimizationFunction(loss, Optimization.AutoZygote())
    elseif autodiff_routine in ["for", "forward"]
        optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
    else
        error("Invalid autodiff routine: $(autodiff_routine)")
    end
    optprob = OptimizationProblem(optf, p_init)
    sol = solve(optprob, Adam(10); maxiters = 1000)

    # Return the parameter estimates
    print("Est\t")
    for k_on in sol
        print(k_on)
        print("\t")
    end
    println()
    
    return sol.u
end

function plot_rn(rn::ReactionSystem,
                 topology::String,
                 p::Vector;
                 plot_type::String="line",
                 t_span::Tuple=(1e-4, 1e5),
                 noise::Float64=0.,
                 num_samples::Int64=100,
                 linesize::Int64=4)::Tuple{Vector, Array, Plots.Plot}

    sample_times, sample_vals, _ = 
        sample_rn(rn, topology, p, t_span, noise=noise, num_samples=num_samples)

    # TODO: THIS CURRENTLY ONLY WORKS FOR FULLY-CONNECTED HETEROTRIMER (N=3) DUE
    # TO THE HARDCODED INDEXING
    if plot_type == "line"
        pl = plot(log10.(sample_times), sample_vals[1, :], lw=linesize, label="Monomer")
        plot!(log10.(sample_times), sample_vals[3, :], lw=linesize, label="Dimer")
        plot!(log10.(sample_times), sample_vals[7, :], lw=linesize, label="Trimer")
    elseif plot_type == "dot"
        pl = scatter(log10.(sample_times), sample_vals[1, :], lw=linesize, label="Monomer")
        scatter!(log10.(sample_times), sample_vals[3, :], lw=linesize, label="Dimer")
        scatter!(log10.(sample_times), sample_vals[7, :], lw=linesize, label="Trimer")
    end
    plot!(xlabel="log10(t)", ylabel="Relative yield")

    return sample_times, sample_vals, pl
end

function get_one_index_per_step(rn::ReactionSystem)::Vector{Int64}
    species_lens = collect(length.(string.(rn.species)))
    
    return [findfirst(len .== species_lens)
            for len in unique(species_lens)]
end