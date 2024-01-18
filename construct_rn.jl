# Code for Julia pkg command line environment setup:
# add Catalyst DifferentialEquations Plots Zygote SciMLSensitivity Optimization OptimizationOptimisers ForwardDiff

using Catalyst
using DifferentialEquations
using Plots
using Zygote
using SciMLSensitivity
using Optimization
using OptimizationOptimisers
#using ForwardDiff

# fc = fully connected
# rt = ring topology

function get_fc_species_list(n::UInt8)::Vector{Vector{String}}
    """
    Given integer `n`, return a length-sorted Vector{Vector{String}} object in which
    Vector{String} objects are all ways of listing a non-empty subset of the `n` 
    monomers in a fully connected reaction network topology, up to permutation.
    Each monomer is represented by a string of the form `Xi`.
    """
    monomers = ["X$i" for i in 1:n]
    species_list = []

    for i in 1:(2 ^ n - 1)
        # `inclusion` uses binary encodings to iterate through non-empty elements
        # of the power set of the set of monomers.
        inclusion = digits(i, base=2, pad=n)
        addend_species = [monomers[i] for i in 1:n if inclusion[i] == 1]
        sort!(addend_species, by=x->parse(Int, x[2:end]))
        species_list = [species_list; [addend_species]]
    end

    sort!(species_list, by=x->length(x))
    return species_list
end

function get_species_string(species_list::Vector{Vector{String}})::String
    """
    Given Vector{Vector{String}} `species_list`, return a String to be evaluated
    as code that lists all species as functions of time for Catalyst.jl's `@species`
    macro. 
    """
    # Convert Vector{String} objects that represent species into Strings
    species_string_list = [join(species) for species in species_list]
    # Prepend macros
    eval_string_line1 = "@variables t; "
    eval_string_line2 = "@species " * join(species_string_list, "(t) ") * "(t)"

    return eval_string_line1 * eval_string_line2
end

function convert_vec_to_string(vec::Vector{String})::String
    """
    Given a Vector object `vec`, return a String representation of it. E.g. given
    a vector containing "X1" and "X2", return String "[X1, X2]".
    """
    return "[" * join(vec, ", ") * "]"
end

function get_fc_rxs_eval_string(species_list::Vector{Vector{String}})::String
    """
    Given a Vector{Vector{String}} `species_list`, return a String to be evaluated
    as code that lists all reactions in a fully connected reaction network topology.
    """
    # In the fully connected topology, the number of species is 2 ^ n - 1.
    num_monomers = Int(log2(length(species_list) + 1))
    # Remove the end product from species_list
    sort!(species_list, by=x->length(x))
    pop!(species_list)

    monomers = [monomer[1] for monomer in species_list[1:num_monomers]]

    rxs_eval_string = "rxs = ["
    # reactant is a vector of bound monomers
    for reactant in species_list
        k_index = length(reactant) * 2 - 1
        # monomer is a singleton vector containing a monomer
        for monomer in monomers
            # Don't record reactions in which either (1) a complex containing itself
            # or (2) the reaction between monomers has already been recorded, just
            # in the other order of reactants.
            if monomer in reactant ||
                (length(reactant) == 1 && Int(monomer[2]) < Int(reactant[1][2]))
                continue
            end

            reactants = [join(reactant), monomer]
            product_string = join(sort([reactant; monomer], by=x->parse(Int, x[2:end])))
            product_string = "[" * product_string * "]"
            reactants_string = convert_vec_to_string(reactants)
            # Add forward direction of the reaction
            rxs_eval_string *= "Reaction(k[$(k_index)], $(reactants_string), $(product_string)), "
            # Add reverse direction of the reaction
            rxs_eval_string *= "Reaction(k[$(k_index+1)], $(product_string), $(reactants_string)), "
        end
    end

    # Prepend macro that specifies the parameters, which we store in the Vector{Float64}
    # `k`
    num_params = (num_monomers - 1) * 2
    eval_string_line1 = "@parameters k[1:$(num_params)]; "
    # Append closing bracket to Vector{Reaction} object
    eval_string_line2 = rxs_eval_string * "]"

    return eval_string_line1 * eval_string_line2
end

function get_rate_constants_from_k_ons(k_ons::Vector,
                                       topology::String;
                                       delta_G_kb_T::Float64=-20., 
                                       C0::Float64=1e4)::Vector{Float64}
    """
    Given binding reaction rates `k_ons`, return a Vector of alternating k_on and
    corresponding k_off values, according to eq. 2 in the preprint.
    """
    rates = []
    for (i, k_on) in enumerate(k_ons)
        if topology == "fc"
            m = i
        elseif topology == "rt"
            m = min(2, i)
        else
            # In the ring topology, 
            error("Topology must be 'fc' or 'rt'")
        end
        rates = [rates; [k_on, k_on * C0 * exp(m * delta_G_kb_T)]]
    end
    
    return rates
end

function get_fc_rn(n::UInt8)::ReactionSystem
    """
    Given integer `n`, return a ReactionSystem object representing a fully connected
    rate growth chemical reaction network.
    """
    species_list = get_fc_species_list(n)

    # Evaluate meta-programmed strings describing the species and reactions involved
    # in the network.
    species_eval_string = get_species_string(species_list)
    eval(Meta.parse(species_eval_string))
    rxs_eval_string = get_fc_rxs_eval_string(species_list)
    eval(Meta.parse(rxs_eval_string))

    @named rn = ReactionSystem(rxs, t)
    return rn
end

function sample_rn(rn::ReactionSystem,
                   topology::String;
                   k_ons::Vector=[10., 10.], 
                   t_start::Float64=1e-4,
                   t_end::Float64=1e4,
                   noise::Float64=0.1)::Tuple{Float64, Vector{Float64}, Array{Float64, 2}, ODEProblem}
    # Start with concentration 1 for all monomers, whose names have 5 characters,
    # e.g. "X1(t)"
    u0 = [length(string(species)) == 5 for species in species(rn)]
    # The time span should start from a small, but non-zero time in order to avoid 
    # taking the log of zero.
    tspan = (t_start, t_end)
    # Obtain a log10-spaced vector of times at which to sample the ODE solution.
    # The log10-scale enables us to observe trapping at early simulation times.
    sample_times = 10 .^ range(log10(tspan[1]), log10(tspan[2]), length=100)

    params = get_rate_constants_from_k_ons(k_ons, topology)
    prob = ODEProblem(rn, u0, tspan, params)
    sol_real = solve(prob, Tsit5(); tstops = sample_times)
    sample_vals = Array(sol_real(sample_times))

    # Add noise to the ODE solution by uniformly sampling a coefficient for the
    # observed value at each time step from the interval [1-`noise`/2, 1+`noise`/2].
    sample_vals .*= (1 .+ noise * rand(Float64, size(sample_vals)) .- (noise / 2))

    return t_end, sample_times, sample_vals, prob
end

function optimise_p(p_init::Vector{Float64}, 
                    topology::String,
                    t_end::Float64, 
                    sample_times::Vector{Float64}, 
                    sample_vals::Array{Float64, 2}, 
                    prob::ODEProblem)::Vector{Float64}
    """
    `topology` should be "fc" for fully connected topology or "rg" for ring topology
    """
    function loss(p, _)
        newtimes = filter(<=(t_end), sample_times)
        params = get_rate_constants_from_k_ons(p, topology)
        newprob = remake(prob; tspan=(1e-4, t_end), p=params)
        sol = Array(solve(newprob, TRBDF2(); saveat = newtimes))
        loss = sum(abs2, sol[:, 7] .- sample_vals[:, 7])
        return loss, sol
    end

    # optimize for the parameters that minimize the loss
    optf = OptimizationFunction(loss, Optimization.AutoZygote())
    optprob = OptimizationProblem(optf, p_init)
    sol = solve(optprob, OptimizationOptimisers.ADAM(10); maxiters = 1000)

    # Return the parameter estimates
    return sol.u
end

# Load reaction network and get sample to be fit to
rn = get_fc_rn(3)
t_end, sample_times, sample_vals, prob = sample_rn(rn, "fc", k_ons=[100, 100], noise=0.)
scatter(log10.(sample_times), sample_vals[7, :])

# Estimate parameters up to `t_end_est`
t_end_est = 1e0
p_init = [10.0, 100.0]
estimate_sol = optimise_p(p_init, "fc", t_end_est, sample_times, sample_vals, prob)     

# Plot the sample and the estimate
scatter(log10.(sample_times), sample_vals[7, :])
_, samples_times_est, sample_vals_est, _ = 
    sample_fc_rn(rn, [estimate_sol[1], 1e-4, estimate_sol[2], 1e-4], t_end_est, 0.0)
plot!(log10.(samples_times_est), sample_vals_est[7, :], lw=4)

# Things to test:
# - Different ODE solvers 
#   - Hypothesis: No significant effect
# - FowardDiff vs. Zygote (reverse-mode) vs. NM (numerical black-box) vs. 
#   gradient descent (numerical first-order) vs. BFGS (numerical second-order)
#   - Hypothesis: Zygote is fastest, but ForwardDiff avoids memory issues at high `n`
# - Different `n` values, e.g. {3, 4, 5, 10}
# - Fully connected vs. ring topologies
# - Optimization using observations of different species