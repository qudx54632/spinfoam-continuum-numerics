"""
Orientation convention for ordered 4-simplices.

A positive oriented simplex is written as its ordered vertex list

    (v1, v2, v3, v4, v5).

A negative oriented simplex is represented by swapping the first two vertices:

    -(v1, v2, v3, v4, v5)  ->  (v2, v1, v3, v4, v5).

No extra sign is multiplied later; the vertex order itself carries the
orientation.
"""

"""Return one fixed representative of `sign * simplex`."""
function simplex_with_orientation(simplex, sign::Integer)
    v = Tuple(simplex)

    if sign == 1
        return v
    elseif sign == -1
        return (v[2], v[1], v[3], v[4], v[5])
    else
        throw(ArgumentError("orientation sign must be +1 or -1"))
    end
end

"""The oriented 4-simplex vertex order used by the local magnetic convention."""
oriented_four_simplex(simplex) = Tuple(simplex)

"""
Five tetrahedra of an ordered 4-simplex, in opposite-vertex order.

For `(v1,v2,v3,v4,v5)`, the output is

    opposite v1, opposite v2, ..., opposite v5.

The tetrahedron keys themselves are sorted, because the same geometric
tetrahedron must have the same dictionary key in neighbouring 4-simplices.
"""
function oriented_simplex_tetrahedra(simplex)
    v = oriented_four_simplex(simplex)

    return (
        tetrahedron_key(v[2], v[3], v[4], v[5]),
        tetrahedron_key(v[1], v[3], v[4], v[5]),
        tetrahedron_key(v[1], v[2], v[4], v[5]),
        tetrahedron_key(v[1], v[2], v[3], v[5]),
        tetrahedron_key(v[1], v[2], v[3], v[4]),
    )
end

"""
Ten triangle keys in the local magnetic spin order.

For the ordered 4-simplex `(v1,v2,v3,v4,v5)`, the spin order is

    (j12, j13, j14, j15, j23, j24, j25, j34, j35, j45),

where `jab` labels the triangle opposite the local pair `(va,vb)`.
The triangle keys are sorted for dictionary lookup.
"""
function oriented_simplex_spin_keys(simplex)
    v = oriented_four_simplex(simplex)

    opposite_triangle(a, b) =
        triangle_key((v[c] for c in 1:5 if c != a && c != b)...)

    return (
        opposite_triangle(1, 2),
        opposite_triangle(1, 3),
        opposite_triangle(1, 4),
        opposite_triangle(1, 5),
        opposite_triangle(2, 3),
        opposite_triangle(2, 4),
        opposite_triangle(2, 5),
        opposite_triangle(3, 4),
        opposite_triangle(3, 5),
        opposite_triangle(4, 5),
    )
end

"""The five tetrahedron keys in local magnetic intertwiner order."""
oriented_simplex_intertwiner_keys(simplex) = oriented_simplex_tetrahedra(simplex)

"""Read the ten local 4-simplex spins from `spin_data`."""
function oriented_simplex_spins(simplex, spin_data::AbstractDict)
    return Tuple(spin_data[key] for key in oriented_simplex_spin_keys(simplex))
end
