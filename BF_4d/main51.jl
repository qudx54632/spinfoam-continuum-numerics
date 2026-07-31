include(joinpath(@__DIR__, "core", "load_bf4d.jl"))
include(joinpath(@__DIR__, "pachner_moves", "one_five_move.jl"))

# 1-5 Pachner move example.
#
# This example uses the orientation-aware recoupled vertex with the local
# 15j amplitude.  The oriented 1-5 geometry is written explicitly in
# `one_five_move.jl`.

use_fast_6j = true
state_sum_function =
    use_fast_6j ? state_sum_amplitude_recoupled_oriented_6j :
                  state_sum_amplitude_recoupled_oriented

boundary_spins = (
    1 // 1,
    1 // 1,
    1 // 2,
    1 // 2,
    1 // 1,
    1 // 2,
    1 // 2,
    1 // 2,
    1 // 2,
    1 // 2,
)

boundary_intertwiners = (
    1 // 1,
    1 // 1,
    0 // 1,
    0 // 1,
    0 // 1,
)

# boundary_spins = fill(3 // 2, 10)
# boundary_intertwiners = fill(0 // 1, 5)

# Fixed internal triangle spins.  The remaining six internal triangle spins
# are summed up to `j_max`.
fixed_internal_spins = (1//2, 0 // 1, 1//2, 0//1)
j_max = 5 // 2

before = one_five_before_amplitude_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    state_sum_function = state_sum_function,
)

after_fixed = one_five_after_sum_fixed_internal_triangles_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    fixed_internal_spins,
    j_max,
    state_sum_function = state_sum_function,
)

fixed_spin_factor = one_five_fixed_internal_triangle_factor(fixed_internal_spins)
boundary_orientation_phase = pachner_boundary_orientation_phase(boundary_spins)
after_fixed_in_before_convention = boundary_orientation_phase * after_fixed

println("1-5 Pachner move")
println("vertex convention        = orientation-aware recoupled, ", use_fast_6j ? "fast 6j local vertex and recoupling" : "3j local vertex")
println("before simplices         = ", (one_five_before_simplex_oriented,))
println("after simplices          = ", one_five_after_simplices_oriented)
println("boundary spins           = ", boundary_spins)
println("boundary intertwiners    = ", boundary_intertwiners)
println("fixed internal triangles = ", one_five_fixed_internal_triangles)
println("fixed internal spins     = ", fixed_internal_spins)
println("free internal triangles  = ", one_five_free_internal_triangles)
println("j_max                    = ", j_max)
println("before amplitude         = ", before)
println("after fixed-spin sum     = ", after_fixed)
println("boundary phase           = ", boundary_orientation_phase)
println("fixed-spin factor        = ", fixed_spin_factor)
println("phase-corrected after / factor = ", after_fixed_in_before_convention / fixed_spin_factor)
println("difference               = ", after_fixed_in_before_convention / fixed_spin_factor - before)
