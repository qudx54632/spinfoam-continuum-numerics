include(joinpath(@__DIR__, "..", "semiclassical", "semiclassical_vertex_tools.jl"))

# Compare the SU(2) BF coherent single-vertex amplitude computed in two ways:
#
#   1. the reference magnetic 3j contraction;
#   2. the fast 6j local vertex plus fast 6j recoupling.
#
# The boundary state is a Livine-Speziale coherent boundary state.  The normal
# data below are generated so that each boundary tetrahedron satisfies closure.

simplex = (1, 2, 3, 4, 5)

# Standard equispin scaling:
#
#     j_ab(lambda) = lambda/2.
#
# Change this integer directly, or run for example
#
#     BF_COHERENT_LAMBDA=3 julia checks/check_su2_bf_coherent_3j_vs_6j_timing.jl
lambda = parse(Int, get(ENV, "BF_COHERENT_LAMBDA", "2"))
j = lambda // 2

spin_data, normal_data = regular_4simplex_boundary_data(simplex, j)
boundary_spins = Tuple(
    spin_data[key]
    for key in oriented_simplex_spin_keys(simplex)
)
closure_error = max_closure_norm(simplex, spin_data, normal_data)

# Warm-up both code paths first.  This avoids mixing Julia compilation time
# into the numerical timing.
coherent_su2_bf_single_vertex_amplitude(
    simplex,
    spin_data,
    normal_data,
)

coherent_su2_bf_single_vertex_amplitude_6j(
    simplex,
    spin_data,
    normal_data,
)

time_3j = @elapsed W_3j = coherent_su2_bf_single_vertex_amplitude(
    simplex,
    spin_data,
    normal_data,
)

time_6j = @elapsed W_6j = coherent_su2_bf_single_vertex_amplitude_6j(
    simplex,
    spin_data,
    normal_data,
)

println("SU(2) BF coherent vertex: 3j reference vs fast 6j")
println("simplex          = ", simplex)
println("lambda           = ", lambda)
println("j_ab             = ", j)
println("spin order       = (j12,j13,j14,j15,j23,j24,j25,j34,j35,j45)")
println("boundary spins   = ", boundary_spins)
println("closure error    = ", closure_error)
println()
println("W_3j             = ", W_3j)
println("W_6j             = ", W_6j)
println("W_6j - W_3j      = ", W_6j - W_3j)
println("|W_6j - W_3j|    = ", abs(W_6j - W_3j))
println()
println("time 3j backend  = ", time_3j, " seconds")
println("time 6j backend  = ", time_6j, " seconds")
println("speedup 3j/6j    = ", time_3j / time_6j)
