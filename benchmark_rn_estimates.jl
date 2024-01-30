using BenchmarkTools

include("construct_rn.jl")
include("estimate_rn_params.jl")

# Kinetic traps will need to be fed in to optimization in a special way... 
# I need a non-autodiff method, ideally both NM and GD/BFGS, as (a) control(s)
function sample_k_ons(rn::ReactionSystem, n::Int64)::Vector
    is_trapped = false
    k_ons = rand(Float64, (n-1)) .* 50. .+ 100.

    while !is_trapped
        k_ons = rand(Float64, (n-1)) .* 50. .+ 100.
        st = log10_range(1e-4, 1e5)
        sample_vals, _ = sample_rn(rn, "fc", k_ons, st)
        kinetic_traps = detect_traps(st, sample_vals[end, :])
        is_trapped = length(kinetic_traps) > 0
    end

    return k_ons
end

# Function to be benchmarked
function optimization_benchmarking_wrapper(n::Int64, ns::Int64, ar::String="reverse", et::String="end_product_yield")
    rn = get_fc_rn(n)

    top = "fc"
    ts = (1e-4, 1e5)
    st = log10_range(ts[1], ts[2], ns)

    b = @benchmarkable optimise_p(p_init, top, st, sv_and_prob[1], sv_and_prob[2], ar, et) setup=
        (st=$st; 
         sv_and_prob = $sample_rn($rn, $top, $sample_k_ons($rn, $n), $st, noise=0., printer=true);
         p_init=rand(Float64, ($n-1)) .* 50. .+ 100.; 
         top=$top; 
         ar=$ar;
         et=$et)
 
    return b
end

for n in [3, 4, 5]
    for mode in ["forward", "reverse"]
        if n != 5
            continue
        end

        print("n\t")
        print(n)
        println()
        print("mode\t")
        print(mode)
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


# Bounds [x∈[1, ]]

# Time optimise_p while varying
# - `p_init` (uniform random sampling from ????)
# - `prob` (object returned by `ODEProblem(rn, u0, t_span, params)`)
#   - by varying 
#     - `rn`
#       - by varying `n` in, e.g., {3, 4, 5, 10}
#     - `params` (uniform random sampling from ...)


# t_span_est = (1e-4, 1e5)
# p_init = [10.0, 10.0]
#estimate_sol = optimise_p(p_init, "fc", t_span_est, sample_times, sample_vals, prob)     


# Throwaway testing wrapper
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

# rn = get_fc_rn(3)
# plot_and_traps(rn, rand(Float64, (2)) .* 100 .+ 25)

# rn = get_fc_rn(4)
# plot_and_traps(rn, rand(Float64, (3)) .* 200 .+ 50)

# rn = get_fc_rn(5)
# plot_and_traps(rn, rand(Float64, (4)) .* 100 .+ 25)

# plot_and_traps(rn, [100, 100])
# plot_and_traps(rn, [200, 200])
# plot_and_traps(rn, [100, 10])
# plot_and_traps(rn, [9075, 100])