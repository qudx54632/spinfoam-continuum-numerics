include("tetrahedron_amplitude.jl")
include("state_sum_utils.jl")
include("one_four_move.jl")
include("one_four_cutoff_plot.jl")

boundary_spins = (1, 1, 1, 1, 1, 1)
j_max = 2
fixed_j04 = 1
max_plot_j_max = 20
cutoffs = 1:0.5:max_plot_j_max

before_tetrahedron = one_four_before_tetrahedron
after_tetrahedra = one_four_after_tetrahedra

before_move_amplitude = one_four_before_amplitude(boundary_spins)
admissible_internal_spins = one_four_admissible_internal_spins(boundary_spins, j_max)
after_move_sum = one_four_after_sum(boundary_spins, j_max)
fixed_j04_sum = one_four_after_sum_fixed_j04(boundary_spins, fixed_j04, j_max)
fixed_j04_physical = one_four_physical_amplitude_fixed_j04(boundary_spins, fixed_j04, j_max)
cutoff_data = one_four_cutoff_comparison_data(boundary_spins, cutoffs)
plot_x = [row.j_max for row in cutoff_data]
after_sum_values = [row.after_sum for row in cutoff_data]
be_prediction_values = [row.be_prediction for row in cutoff_data]

println("1-4 move before tetrahedron = ", before_tetrahedron)
println("1-4 move after tetrahedra = ", after_tetrahedra)
println("before move amplitude     = ", Float64(before_move_amplitude))
println("j_max                     = ", j_max)
println("number of admissible terms = ", length(admissible_internal_spins))
println("after cutoff sum amplitude = ", Float64(after_move_sum))
println("fixed j04                 = ", fixed_j04)
println("fixed j04 sum             = ", Float64(fixed_j04_sum))
println("fixed j04 sum / d_j04^2   = ", Float64(fixed_j04_physical))
println("fixed j04 test difference = ", Float64(fixed_j04_physical - before_move_amplitude))
println("plot cutoff range         = ", first(cutoffs), ":", step(cutoffs), ":", last(cutoffs))
println("max plot j_max            = ", max_plot_j_max)
println("last plot ratio           = ", after_sum_values[end] / be_prediction_values[end])


using Plots
plot(plot_x, after_sum_values;
    xscale = :log10,
    yscale = :log10,
    marker = :circle,
    label = "after cutoff sum",
)
plot!(plot_x, be_prediction_values;
    marker = :square,
    linestyle = :dash,
    label = "BE factor * before",
)
