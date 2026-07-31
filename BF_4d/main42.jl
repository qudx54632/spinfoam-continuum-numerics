include(joinpath(@__DIR__, "core", "load_bf4d.jl"))
include(joinpath(@__DIR__, "pachner_moves", "two_four_move.jl"))

# 2-4 Pachner move example.
#
# This example uses the orientation-aware recoupled vertex with the local
# 15j amplitude.  The oriented 2-4 geometry is written explicitly in
# `two_four_move.jl`.

use_fast_6j = true
state_sum_function =
    use_fast_6j ? state_sum_amplitude_recoupled_oriented_6j :
                  state_sum_amplitude_recoupled_oriented

boundary_spins = (
    2 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
    1 // 1,
)

boundary_intertwiners = (
    1 // 1,
    0 // 1,
    0 // 1,
    0 // 1,
    1 // 1,
    0 // 1,
    0 // 1,
    0 // 1,
)

fixed_position = 1
fixed_spin = 2 // 1
j_max = 4

before = two_four_before_sum_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    state_sum_function = state_sum_function,
)

after_fixed = two_four_after_sum_fixed_internal_triangle_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    fixed_position,
    fixed_spin,
    j_max,
    state_sum_function = state_sum_function,
)

fixed_spin_factor = dimension(fixed_spin)^2
boundary_orientation_phase = pachner_boundary_orientation_phase(boundary_spins)
after_fixed_in_before_convention = boundary_orientation_phase * after_fixed

println("2-4 Pachner move")
println("vertex convention        = orientation-aware recoupled, ", use_fast_6j ? "fast 6j local vertex and recoupling" : "3j local vertex")
println("before simplices         = ", two_four_before_simplices_oriented)
println("after simplices          = ", two_four_after_simplices_oriented)
println("boundary spins           = ", Float64.(boundary_spins))
println("boundary intertwiners    = ", Float64.(boundary_intertwiners))
# println("internal triangles       = ", two_four_after_internal_triangles)
println("fixed internal triangle  = ", two_four_after_internal_triangles[fixed_position])
println("fixed spin               = ", fixed_spin)
println("j_max                    = ", j_max)
println("before amplitude         = ", before)
println("after fixed-spin sum     = ", after_fixed)
println("boundary phase           = ", boundary_orientation_phase)
println("fixed-spin factor        = ", fixed_spin_factor)
println("phase-corrected after / factor = ", after_fixed_in_before_convention / fixed_spin_factor)
println("difference               = ", after_fixed_in_before_convention / fixed_spin_factor - before)
