using Catalyst
using DifferentialEquations
using Optimization
using OptimizationOptimisers
using ForwardDiff
using Zygote
using SciMLSensitivity
using OptimizationOptimJL

include("construct_rn.jl")

"""
    log10_range(start, ending[, length])

Like the built-in `range()` function, but returns values that are equidistant in
log10-space.
"""
function log10_range(start, ending, length::Int64=100)::Vector
    return 10 .^ range(log10(start), log10(ending), length=length)
end

"""
    sample_rn(rn, topology, k_ons, sample_times[, noise, printer])

Solve the reaction network `rn` parameterized by `topology` and `k_ons` and return 
a "sample" of its concentrations at times listed in `sample_times`.
"""
function sample_rn(rn::ReactionSystem,
                   topology::String,
                   k_ons::Vector,
                   sample_times::Vector;
                   noise::Float64=0.,
                   printer::Bool=false)

    # Print output true k_on values. Used to get output from benchmarking in 
    # benchmark_run_estimates.jl
    if printer
        print("True\t")
        for k_on in k_ons
            print(k_on, "\t")
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

"""
    kinetic_trap(y, t_start, t_end)

Custom type used to record the y-value and horizontal bounds of a kinetic trap.
"""
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
`threshold=0.05`. Flat reginos of `u` at the beginning and end of the provided data
will not be counted as traps, as they represent initial and terminal equilibrium,
not intermediate pseudo-equilibriated states.
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

        # Prevent counting the beginning of a reaction, which will likely be flat
        # as a trap
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
                    prob::ODEProblem;
                    autodiff_routine::String="reverse",
                    experiment_type::String="end_product_yield",
                    printer::Bool=true)::Vector{Float64}
    # `topology` should be "fc" for fully connected topology or "rg" for ring topology

    if experiment_type == "kinetic_traps"
        # Ensure 100 sample times for trap detection
        sample_times = log10_range(sample_times[1], sample_times[end])
        # sample_vals should be noise-free
        kinetic_traps = detect_traps(sample_times, sample_vals[end, :])
    end

    # This weird generic typing is necessitated by the explicit use of dual numbers
    # by forward-mode autodiff, as implemented in `Optimization.AutoForwardDiff()`.
    function loss(p::AbstractVector{T}, _) where T
        params = get_rate_constants_from_k_ons(p, topology)

        # params = T[params[1],
        #            params[2],
        #            params[3],
        #            params[4],
        #            params[5],
        #            params[6],
        #            params[7],
        #            params[8]]

        newprob = remake(prob; 
                         p=params, 
                         tspan=T[sample_times[1], sample_times[end]])
        # TRBDF2 is a stiff solver that is warranted by the log-time-scale of kinetic
        # kinetic traps, there are situations in which Tsit5 can be used instead.
        sol = Array(solve(newprob, TRBDF2(); saveat=sample_times, maxiters=1e7))

        if experiment_type == "end_product_yield"
            loss = sum((sol[end, :] .- sample_vals[end, :]) .^ 2)
        elseif experiment_type == "kinetic_traps"
            # Calculate mean-squared error between the theoretical concentration 
            # trajectory and the concentrations of observed kinetic traps. Importantly,
            # this loss function does not take into account any information about 
            # the concentrations of the species at times at which kinetic traps
            # were not observed.
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
    if autodiff_routine in ["r", "rev", "reverse"]
        optf = OptimizationFunction(loss, Optimization.AutoZygote())
    elseif autodiff_routine in ["f", "for", "forward"]
        optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
    else
        error("Invalid autodiff routine: $(autodiff_routine)")
    end

    optprob = OptimizationProblem(optf, p_init)
    # The argument passed into `Adam()` is its learning rate.
    sol = solve(optprob, Adam(10); maxiters=1000)

    # Return the parameter estimates
    if printer
        print("Est\t")
        for k_on in sol
            print(k_on, "\t")
        end
        println()
    end

    return sol.u
end

"""
    get_one_index_per_step(rn)

Returns a vector of the indices of the first (lexicographically) species at each
level of assembly (e.g. monomer, dimer, trimer...) in the reaction network `rn`.
This function is not currently used, but could be used to run an optimization routine
that observes one species at each level of assembly. The motivation for this is
that, in the rate growth model, all species at a given level of assembly have the
same concentration at all times. 
"""
function get_one_index_per_step(rn::ReactionSystem)::Vector{Int64}
    species_lens = collect(length.(string.(rn.species)))
    
    return [findfirst(len .== species_lens)
            for len in unique(species_lens)]
end