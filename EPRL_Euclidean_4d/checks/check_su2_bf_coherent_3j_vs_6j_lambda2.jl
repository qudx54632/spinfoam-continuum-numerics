include(joinpath(@__DIR__, "..", "semiclassical", "semiclassical_vertex_tools.jl"))

simplex = (1, 2, 3, 4, 5)

# Boundary spins in the local 4-simplex order
#
#     (j12,j13,j14,j15,j23,j24,j25,j34,j35,j45),
#
# where jab labels the triangle opposite the local vertex pair (a,b).
# Change these ten numbers directly.
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

spin_data = Dict(
    key => spin
    for (key, spin) in zip(oriented_simplex_spin_keys(simplex), boundary_spins)
)

# Construct one closed set of coherent normals for each boundary tetrahedron.
# This guarantees tetrahedron-wise closure for non-equal spins.
normal_data = closed_boundary_normal_data(simplex, spin_data)

closure_error = max_closure_norm(simplex, spin_data, normal_data)

println("SU(2) BF coherent vertex: old 3j backend vs fast 6j backend")
println("simplex        = ", simplex)
println("spin order     = (j12,j13,j14,j15,j23,j24,j25,j34,j35,j45)")
println("boundary spins = ", boundary_spins)
println("closure error  = ", closure_error)
println()

W_3j = coherent_su2_bf_single_vertex_amplitude(simplex, spin_data, normal_data)
W_6j = coherent_su2_bf_single_vertex_amplitude_6j(simplex, spin_data, normal_data)

println("W old 3j       = ", W_3j)
println("W fast 6j      = ", W_6j)
println("difference     = ", W_6j - W_3j)
println("|difference|   = ", abs(W_6j - W_3j))
