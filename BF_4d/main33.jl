include(joinpath(@__DIR__, "core", "load_bf4d.jl"))
include(joinpath(@__DIR__, "pachner_moves", "three_three_move.jl"))

# 3-3 Pachner move example.
#
# This example uses the orientation-aware recoupled vertex with the local
# 15j amplitude.  The oriented 3-3 geometry is written explicitly in
# `three_three_move.jl`.

use_fast_6j = true
state_sum_function =
    use_fast_6j ? state_sum_amplitude_recoupled_oriented_6j :
                  state_sum_amplitude_recoupled_oriented

boundary_spins = (
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
    1 // 1,
    1 // 1,
    2 // 1,
)

boundary_intertwiners = (
    0 // 1,
    0 // 1,
    0 // 1,
    1 // 1,
    0 // 1,
    1 // 1,
    1 // 1,
    0 // 1,
    1 // 1,
)

j_max = 4

before = three_three_before_sum_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    j_max,
    state_sum_function = state_sum_function,
)

after = three_three_after_sum_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    j_max,
    state_sum_function = state_sum_function,
)

boundary_orientation_phase = pachner_boundary_orientation_phase(boundary_spins)
after_in_before_convention = boundary_orientation_phase * after

println("3-3 Pachner move")
println("vertex convention       = orientation-aware recoupled, ", use_fast_6j ? "fast 6j local vertex and recoupling" : "3j local vertex")
println("before simplices        = ", three_three_before_simplices_oriented)
println("after simplices         = ", three_three_after_simplices_oriented)
println("boundary spins          = ", Float64.(boundary_spins))
println("boundary intertwiners   = ", Float64.(boundary_intertwiners))
# println("before internal triangle = ", three_three_before_internal_triangles[1])
# println("after internal triangle  = ", three_three_after_internal_triangles[1])
println("j_max                   = ", j_max)
println("before amplitude        = ", before)
println("after amplitude         = ", after)
println("boundary phase          = ", boundary_orientation_phase)
println("phase-corrected after   = ", after_in_before_convention)
println("difference              = ", after_in_before_convention - before)
