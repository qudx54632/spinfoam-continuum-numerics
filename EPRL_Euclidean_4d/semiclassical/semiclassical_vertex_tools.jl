include(joinpath(@__DIR__, "..", "load_eprl.jl"))

using LinearAlgebra

"""
Single-vertex coherent-boundary semiclassical tools.

References:
    - Livine-Speziale coherent intertwiners: arXiv:0705.0674.
    - EPRL/FK coherent vertex asymptotics and Regge action:
      Barrett et al., arXiv:0902.1170.
"""

"""
Regular 4-simplex normals in R^5.

The five unit vectors satisfy

    sum_a N_a = 0,              N_a . N_b = -1/4  for a != b.

Here `N_a` labels the outward normal to the tetrahedron opposite vertex `a`.
"""
function regular_4simplex_normals_4d()
    normals = Vector{Float64}[]
    center = fill(1 / 5, 5)

    for a in 1:5
        e = zeros(5)
        e[a] = 1
        v = e - center
        push!(normals, v / norm(v))
    end

    return normals
end

function hyperplane_frame(Na)
    one = fill(1 / sqrt(5), 5)
    frame = Vector{Float64}[]

    for k in 1:5
        v = zeros(5)
        v[k] = 1

        v -= dot(v, one) * one
        v -= dot(v, Na) * Na

        for q in frame
            v -= dot(v, q) * q
        end

        if norm(v) > 1e-12
            push!(frame, v / norm(v))
        end

        length(frame) == 3 && return frame
    end

    throw(ArgumentError("could not construct a 3D frame"))
end

function local_coordinates(v, frame)
    coords = Tuple(dot(v, q) for q in frame)
    length_coords = sqrt(sum(x^2 for x in coords))

    return Tuple(x / length_coords for x in coords)
end

"""
Regular 4-simplex boundary data.

For each boundary tetrahedron tau and triangle f inside tau, the normal data
are keyed by `(tau,f)`.  This is necessary because the same triangle has
different normals when seen from its two neighbouring tetrahedra.
"""
function regular_4simplex_boundary_data(simplex, j)
    vertices = Tuple(simplex)
    N = regular_4simplex_normals_4d()
    frames = [hyperplane_frame(N[a]) for a in 1:5]

    spin_data = Dict(
        triangle => j
        for triangle in oriented_simplex_spin_keys(simplex)
    )

    normal_data = Dict()

    for a in 1:5
        tau = tetrahedron_key((vertices[c] for c in 1:5 if c != a)...)

        for b in 1:5
            a == b && continue

            f = triangle_key((vertices[c] for c in 1:5 if c != a && c != b)...)
            projected = N[b] - dot(N[b], N[a]) * N[a]
            n = local_coordinates(projected / norm(projected), frames[a])
            normal_data[(tau, f)] = n
        end
    end

    return spin_data, normal_data
end

function max_closure_norm(simplex, spin_data, normal_data)
    max_norm = 0.0

    for tau in oriented_simplex_tetrahedra(simplex)
        spin_keys = global_intertwiner_spin_keys(tau)
        spins = Tuple(spin_data[key] for key in spin_keys)
        normals = Tuple(
            tetrahedron_triangle_normal(normal_data, tau, key)
            for key in spin_keys
        )
        closure = closure_vector(spins, normals)
        max_norm = max(max_norm, norm(collect(closure)))
    end

    return max_norm
end

"""Regular 4-simplex Regge action for equal spins."""
function regular_4simplex_regge_action(j, gamma; dihedral = :interior)
    theta =
        if dihedral == :interior
            acos(1 / 4)
        elseif dihedral == :exterior
            acos(-1 / 4)
        else
            throw(ArgumentError("dihedral must be :interior or :exterior"))
        end

    return Float64(gamma) * 10 * Float64(j) * theta
end

"""Regular 4-simplex SU(2) BF Regge action for equal spins."""
function regular_4simplex_su2_bf_regge_action(j; dihedral = :interior)
    theta =
        if dihedral == :interior
            acos(1 / 4)
        elseif dihedral == :exterior
            acos(-1 / 4)
        else
            throw(ArgumentError("dihedral must be :interior or :exterior"))
        end

    return 10 * Float64(j) * theta
end

function print_semiclassical_vertex_check(
    simplex,
    gamma,
    base_spin,
    lambda_values;
    boundary_data = regular_4simplex_boundary_data,
)
    println("Euclidean EPRL coherent single-vertex semiclassical check")
    println("simplex                = ", simplex)
    println("gamma                  = ", gamma)
    println("base spin              = ", base_spin)
    println()
    println("lambda  j       |W|                 arg(W)       S_int mod 2pi   S_ext mod 2pi   cos(S_int)")

    for lambda in lambda_values
        j = lambda * base_spin
        spin_data, normal_data = boundary_data(simplex, j)

        closure_error = max_closure_norm(simplex, spin_data, normal_data)
        W = coherent_single_vertex_amplitude(simplex, spin_data, normal_data, gamma)

        S_int = regular_4simplex_regge_action(j, gamma; dihedral = :interior)
        S_ext = regular_4simplex_regge_action(j, gamma; dihedral = :exterior)

        println(
            lambda, "       ",
            Float64(j), "     ",
            abs(W), "    ",
            angle(W), "    ",
            mod(S_int, 2pi), "    ",
            mod(S_ext, 2pi), "    ",
            cos(S_int),
        )
        println("        closure error = ", closure_error)
        println("        W             = ", W)
    end

    println()
    println("Expected large-spin structure:")
    println("  W(lambda) ~ lambda^(-12) [N_+ exp(i S_Regge) + N_- exp(-i S_Regge)]")
    println("For a fixed-spin coherent boundary state this is usually a cosine-like")
    println("Regge oscillation, not a single isolated exp(i S_Regge) branch.")
end
