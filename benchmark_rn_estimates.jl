using BenchmarkTools

include("construct_rn.jl")
include("estimate_rn_params.jl")

"""
This file contains functions for constructing, simulating from, and estimating parameters
of reaction networks by referencing `construct_rn.jl` and `estimate_rn_params.jl`.
"""

"""
    sample_k_ons(rn, n[, time_span, topology])

Given a reaction network `rn` and the number of distinct monomers `n` that it contains,
return a Vector of `n-1` random k_on values chosen such that kinetic trapping of
the end product is guaranteed.
TODO: Actually use theory to more efficiently choose parameters that result in
kinetic trapping, instead of blindly iterating until parameters from within a small
range are found.
"""
function sample_k_ons(rn::ReactionSystem, 
                      n::Int64; 
                      time_span::Tuple=(1e-4, 1e5),
                      topology="fc")::Vector
    is_trapped = false

    sample_times = log10_range(time_span[1], time_span[2])

    rand_k_ons() = rand(Float64, (n-1)) .* 50. .+ 100.
    k_ons = rand_k_ons()

    # Repeatedly generate new `k_ons` until choices that result in kinetic trapping
    # are found. This can almost certainly be accelerated by applying theory of
    # kinetic trapping.
    while !is_trapped
        k_ons = rand_k_ons()
        sample_vals, _ = sample_rn(rn, topology, k_ons, sample_times)
        kinetic_traps = detect_traps(sample_times, sample_vals[end, :])
        is_trapped = length(kinetic_traps) > 0
    end

    return k_ons
end

"""
    optimization_benchmarking_wrapper(n, ns[, ar, et])

Running this function performs benchmarking, with output printed to stdout.

"""
function optimization_benchmarking_wrapper(n::Int64, ns::Int64, ar::String="reverse", et::String="end_product_yield")
    rn = get_fc_rn(n)

    top = "fc"
    ts = (1e-4, 1e5)
    st = log10_range(ts[1], ts[2], ns)

    b = @benchmarkable optimise_p(p_init, top, st, sv_and_prob[1], sv_and_prob[2], autodiff_routine=ar, experiment_type=et) setup=
        (st=$st; 
         sv_and_prob = $sample_rn($rn, $top, $sample_k_ons($rn, $n), $st, noise=0., printer=true);
         p_init=rand(Float64, ($n-1)) .* 50. .+ 100.; 
         top=$top; 
         ar=$ar;
         et=$et)
 
    return b
end

for n in [3, 4, 5]
    for mode in ["reverse", "forward"]
        println("n\t", n)
        println("mode\t", mode)
        println()

        b = optimization_benchmarking_wrapper(n, 100, mode, "kinetic_traps")
        results = run(b, seconds=150)

        for time in results.times
            println(time)
        end
        print(results.times)

        println()

        break
    end
end

# Throwaway testing wrapper for `detect_traps()`, should be deleted soon
function plot_and_traps(rn, params)
    println(params)
    sample_times, sample_vals, pl = plot_rn(rn, "fc", params)
    kinetic_traps = detect_traps(sample_times, sample_vals[7, :])

    print(kinetic_traps)

    for kt in kinetic_traps
        plot!(log10.([kt.t_start, kt.t_end]), 
              [kt.y, kt.y], color=:red, linewidth=4)
    end    
    display(pl)
end