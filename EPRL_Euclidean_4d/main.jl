include(joinpath(@__DIR__, "load_eprl.jl"))

# 1-5 Pachner move with the same boundary data before and after refinement.

gamma = 0 // 1

boundary_spins = (
    2 // 1, 1 // 1, 2 // 1, 1 // 1, 1 // 1,
    2 // 1, 1 // 1, 1 // 1, 1 // 1, 1 // 1,
)

boundary_intertwiners = (
    1 // 1, 1 // 1, 0 // 1, 0 // 1, 0 // 1,
)

# Fixed internal triangle spins:
#
#     (012), (013), (014), (015).
#
# This mixed choice avoids the very sharp single-shell jump caused by fixing
# all four internal spins to zero.
fixed_internal_spins = (1 // 1, 1 // 1, 0 // 1, 0 // 1)

# Cutoffs 0, 1/2, 1, ..., 5.
cutoffs = half_integer_spins(5 // 1)

before = one_five_before_amplitude(
    boundary_spins,
    boundary_intertwiners,
    gamma,
)



println("Euclidean EPRL 1-5 Pachner move: fixed internal spins")
println("gamma                  = ", gamma)
println("before simplices       = ", (one_five_before_simplex,))
println("after simplices        = ", one_five_after_simplices)
println("boundary spins         = ", Float64.(boundary_spins))
println("boundary intertwiners  = ", Float64.(boundary_intertwiners))
println("fixed triangles        = ", one_five_fixed_internal_triangles)
println("fixed internal spins   = ", Float64.(fixed_internal_spins))
println("free triangles         = ", one_five_free_internal_triangles)
println("cutoffs                = ", Float64.(cutoffs))
println("allowed spin values    = ", one_five_internal_spin_values(maximum(cutoffs), gamma))
println("before amplitude       = ", before)
