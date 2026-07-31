"""Boundary normal data for coherent intertwiners."""

function pair_vectors_with_resultant(a, b, R, transverse_axis)
    a = Float64(a)
    b = Float64(b)
    R = Float64(R)

    if R == 0.0
        a == b || throw(ArgumentError("zero resultant requires equal lengths"))

        if transverse_axis == :y
            return ((0.0, a, 0.0), (0.0, -a, 0.0))
        elseif transverse_axis == :z
            return ((0.0, 0.0, a), (0.0, 0.0, -a))
        else
            throw(ArgumentError("transverse_axis must be :y or :z"))
        end
    end

    x = (a^2 + R^2 - b^2) / (2R)
    transverse_squared = max(a^2 - x^2, 0.0)
    transverse = sqrt(transverse_squared)

    if transverse_axis == :y
        return ((x, transverse, 0.0), (R - x, -transverse, 0.0))
    elseif transverse_axis == :z
        return ((x, 0.0, transverse), (R - x, 0.0, -transverse))
    else
        throw(ArgumentError("transverse_axis must be :y or :z"))
    end
end

"""
Construct four unit normals for given face areas.

The output normals satisfy

    sum_a areas[a] * normals[a] = 0.

This constructs a closed vector tetrahedron.  It is useful for coherent
boundary-state tests with non-equal spins.  By itself it does not impose
shape matching between different boundary tetrahedra.
"""
function closed_tetrahedron_normals_from_areas(areas)
    length(areas) == 4 || throw(ArgumentError("a tetrahedron has four face areas"))
    all(A -> A > 0, areas) ||
        throw(ArgumentError("closed normal construction currently assumes positive face areas"))

    A = Float64.(areas)
    largest = maximum(A)
    largest <= sum(A) - largest + 1e-12 ||
        throw(ArgumentError("face areas do not satisfy the closure polygon inequality"))

    a, b, c, d = A

    lower = max(abs(a - b), abs(c - d))
    upper = min(a + b, c + d)
    lower <= upper + 1e-12 ||
        throw(ArgumentError("could not choose a common resultant for the four areas"))

    R = (lower + upper) / 2

    if R < 1e-14
        v1, v2 = pair_vectors_with_resultant(a, b, 0.0, :y)
        v3, v4 = pair_vectors_with_resultant(c, d, 0.0, :z)
    else
        v1, v2 = pair_vectors_with_resultant(a, b, R, :y)
        w3, w4 = pair_vectors_with_resultant(c, d, R, :z)
        v3 = (-w3[1], -w3[2], -w3[3])
        v4 = (-w4[1], -w4[2], -w4[3])
    end

    vectors = (v1, v2, v3, v4)

    return Tuple(
        Tuple(vectors[a][mu] / A[a] for mu in 1:3)
        for a in 1:4
    )
end

"""
Construct normal data for all boundary tetrahedra of one 4-simplex.

For each tetrahedron the four normals are constructed independently from the
four spins in the fixed global intertwiner order.
"""
function closed_boundary_normal_data(simplex, spin_data::AbstractDict)
    normal_data = Dict()

    for tetrahedron in oriented_simplex_tetrahedra(simplex)
        spin_keys = global_intertwiner_spin_keys(tetrahedron)
        areas = Tuple(spin_data[key] for key in spin_keys)
        normals = closed_tetrahedron_normals_from_areas(areas)

        for (triangle, normal) in zip(spin_keys, normals)
            normal_data[(tetrahedron, triangle)] = normal
        end
    end

    return normal_data
end
