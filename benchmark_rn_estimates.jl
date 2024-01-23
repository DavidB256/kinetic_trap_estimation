include("construct_rn.jl")
include("estimate_rn_params.jl")

function main()
    # Time optimise_p while varying
    # - `p_init` (uniform random sampling)
    # - `prob` (object returned by `ODEProblem(rn, u0, t_span, params)`)
    #   - by varying 
    #     - `rn`
    #       - by varying `n` in, e.g., {3, 4, 5, 10}
    #     - `params` (uniform random sampling from ...)

    # Load reaction network and get sample to be fit to
    rn = get_fc_rn(3)

    plot_rn(rn, "fc", [100, 100])
    plot_rn(rn, "fc", [200, 200])
    plot_rn(rn, "fc", [100, 10])
    plot_rn(rn, "fc", [10, 100])

    t_span_est = (1e-4, 1e5)
    p_init = [10.0, 10.0]
    #estimate_sol = optimise_p(p_init, "fc", t_span_est, sample_times, sample_vals, prob)     
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end