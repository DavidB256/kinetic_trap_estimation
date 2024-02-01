using Plots

include("construct_rn.jl")

"""
This file contains a script for creating figures that were used in the rotation
talk to illustrate kinetic trapping during heterotrimer assembly and types of 
experimental data.
"""

# Load fully connected heterotrimer reaction network
rn = get_fc_rn(3)

# Set up for plotting
gr(legendfontsize=10)
linesize=4
dotsize=7
xlims_set = (-4.5, 5.5)

# Plot ODE solutions
sample_times, sample_vals, _ = 
    sample_rn(rn, "fc", [100, 100], (1e-4, 1e5), noise=0., num_samples=100)
p1 = plot(log10.(sample_times), sample_vals[1, :], lw=linesize, label="Monomer",
          xlimits=xlims_set, dpi=1000)
plot!(log10.(sample_times), sample_vals[3, :], lw=linesize, label="Dimer")
plot!(log10.(sample_times), sample_vals[7, :], lw=linesize, label="Trimer")
plot!(xlabel="log10(t) (a.u.)", ylabel="Relative yield")
savefig("figs/trimer_fig_pt1.png")

# Plot ODE solution for the end product trimer, but only sampled at 10 time points
sample_times, sample_vals, _ = 
    sample_rn(rn, "fc", [100, 100], (1e-4, 1e5), noise=0., num_samples=10)
p2 = scatter(log10.(sample_times), sample_vals[7, :], label="Trimer",
             xlimits=xlims_set, dpi=1000, ms=dotsize, color=3)
plot!(xlabel="log10(t) (a.u.)", ylabel="Relative yield")
savefig("figs/trimer_fig_pt2.png")

# Create same plot as p2, but with vertical noise added to the points
sample_times, sample_vals, _ = 
    sample_rn(rn, "fc", [100, 100], (1e-4, 1e5), noise=0.2, num_samples=10)
p3 = scatter(log10.(sample_times), sample_vals[7, :], label="Trimer",
             xlimits=xlims_set, dpi=1000, ms=dotsize, color=3)
plot!(xlabel="log10(t) (a.u.)", ylabel="Relative yield")
savefig("figs/trimer_fig_pt3.png")

# Plot only the location and value of the kinetic trap in the end product trimer
p4 = plot([-0.45, 1.05], [0.7, 0.7], lw=linesize, label="Trimer", color=3,
          xlimits=xlims_set, dpi=1000, ylimits=(0., 1.))
plot!(xlabel="log10(t) (a.u.)", ylabel="Relative yield")
savefig("figs/trimer_fig_pt4.png")

# Put all plots into one figure as a grid of 2x2 panels
l = @layout [a b; c d]
plot(p1, p2, p3, p4, layout=l)
savefig("figs/trimer_fig_composite.png")
