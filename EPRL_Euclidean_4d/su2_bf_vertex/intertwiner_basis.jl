"""Fixed global four-valent intertwiner basis for glued tetrahedra."""

"""
Return the four triangle keys of a tetrahedron in one fixed global pairing.

For a sorted tetrahedron `(a,b,c,d)`, the convention is

    ((b,c,d), (a,c,d)) -> i -> ((a,b,d), (a,b,c)).

This is the single intertwiner basis used for all glued tetrahedra.
"""
function global_intertwiner_spin_keys(tetrahedron)
    v = tetrahedron_key(tetrahedron...)

    return (
        triangle_key(v[2], v[3], v[4]),
        triangle_key(v[1], v[3], v[4]),
        triangle_key(v[1], v[2], v[4]),
        triangle_key(v[1], v[2], v[3]),
    )
end

"""Read the four spins entering the fixed global tetrahedron pairing."""
function global_intertwiner_spins(tetrahedron, spin_data::AbstractDict)
    return Tuple(spin_data[key] for key in global_intertwiner_spin_keys(tetrahedron))
end

"""Allowed global intertwiner labels for one tetrahedron."""
function allowed_global_intertwiners(tetrahedron, spin_data::AbstractDict)
    spins = global_intertwiner_spins(tetrahedron, spin_data)
    j1, j2, j3, j4 = spins
    i_max = min(j1 + j2, j3 + j4)

    return [
        i
        for i in half_integer_spins(i_max)
        if admissible_4valent(spins, i)
    ]
end
