"""Return an unoriented edge in a canonical vertex order."""
function edge_key(a::Integer, b::Integer)
    a == b && throw(ArgumentError("an edge needs two different vertices"))
    return a < b ? (a, b) : (b, a)
end

"""Return all half-integer spins from `0` through `j_max`."""
function half_integer_spins(j_max::Real)
    j_max >= 0 || throw(ArgumentError("j_max must be nonnegative"))
    isinteger(2 * j_max) || throw(ArgumentError("j_max must be a half-integer"))
    return [n // 2 for n in 0:Int(2 * j_max)]
end

"""Test the SU(2) admissibility conditions for a triangle of spins."""
function admissible_triangle(spins::Tuple{Real,Real,Real})
    j1, j2, j3 = spins

    return all(j -> j >= 0 && isinteger(2 * j), spins) &&
           j1 + j2 >= j3 &&
           j1 + j3 >= j2 &&
           j2 + j3 >= j1 &&
           isinteger(j1 + j2 + j3)
end

"""Return the six unoriented edges of an ordered tetrahedron."""
function tetrahedron_edges(tetrahedron::Tuple{Integer,Integer,Integer,Integer})
    v1, v2, v3, v4 = tetrahedron

    return (
        edge_key(v1, v2),
        edge_key(v1, v3),
        edge_key(v1, v4),
        edge_key(v2, v3),
        edge_key(v2, v4),
        edge_key(v3, v4),
    )
end

"""Read the three spins around a triangle from `spin_data`."""
function triangle_spins(
    triangle::Tuple{Integer,Integer,Integer},
    spin_data::AbstractDict,
)
    v1, v2, v3 = triangle
    return (
        spin_data[edge_key(v1, v2)],
        spin_data[edge_key(v1, v3)],
        spin_data[edge_key(v2, v3)],
    )
end

"""
Return the two rows of the 6j symbol associated with an ordered tetrahedron.

For `(v1, v2, v3, v4)`, the convention is

    { j12  j13  j23 }
    { j34  j24  j14 }.
"""
function tetrahedron_spin_rows(
    tetrahedron::Tuple{Integer,Integer,Integer,Integer},
    spin_data::AbstractDict,
)
    v1, v2, v3, v4 = tetrahedron

    top = (
        spin_data[edge_key(v1, v2)],
        spin_data[edge_key(v1, v3)],
        spin_data[edge_key(v2, v3)],
    )
    bottom = (
        spin_data[edge_key(v3, v4)],
        spin_data[edge_key(v2, v4)],
        spin_data[edge_key(v1, v4)],
    )

    return top, bottom
end

"""Test whether all four triangular faces of a tetrahedron are admissible."""
function admissible_tetrahedron(
    tetrahedron::Tuple{Integer,Integer,Integer,Integer},
    spin_data::AbstractDict,
)
    v1, v2, v3, v4 = tetrahedron
    faces = (
        (v1, v2, v3),
        (v1, v2, v4),
        (v1, v3, v4),
        (v2, v3, v4),
    )

    return all(face -> admissible_triangle(triangle_spins(face, spin_data)), faces)
end

"""Evaluate the 6j amplitude of a tetrahedron using `spin_data`."""
function tetrahedron_amplitude(
    tetrahedron::Tuple{Integer,Integer,Integer,Integer},
    spin_data::AbstractDict,
)
    top, bottom = tetrahedron_spin_rows(tetrahedron, spin_data)
    return tetrahedron_amplitude(top, bottom)
end

"""Ponzano-Regge edge dimension `d_j = 2j + 1`."""
edge_amplitude(j::Real) = 2 * j + 1

"""Multiply the 6j amplitudes of a collection of tetrahedra."""
function tetrahedra_amplitude(tetrahedra, spin_data::AbstractDict)
    return prod(tetrahedron_amplitude(tetrahedron, spin_data) for tetrahedron in tetrahedra)
end

"""Multiply `2j + 1` over a specified collection of edges."""
function edge_amplitude_product(edges, spin_data::AbstractDict)
    return prod(edge_amplitude(spin_data[edge]) for edge in edges)
end

"""Count how many tetrahedra contain a given edge."""
function edge_incidence(edge, tetrahedra)
    return count(tetrahedron -> edge in tetrahedron_edges(tetrahedron), tetrahedra)
end

"""Ponzano-Regge sign `(-1)^chi` with `chi = sum((n_e - 2) j_e) + chi0`."""
function state_sum_sign(tetrahedra, spin_data::AbstractDict, internal_edges; chi0 = 0)
    chi = chi0 + sum(
        (edge_incidence(edge, tetrahedra) - 2) * spin_data[edge]
        for edge in internal_edges
    )
    isinteger(chi) || throw(ArgumentError("the state-sum sign exponent chi must be an integer"))

    return (-1)^Int(chi)
end

"""Evaluate the fixed-spin Ponzano-Regge summand of a simplicial complex."""
function state_sum_amplitude(tetrahedra, spin_data::AbstractDict, internal_edges; chi0 = 0)
    return state_sum_sign(tetrahedra, spin_data, internal_edges; chi0 = chi0) *
           edge_amplitude_product(internal_edges, spin_data) *
           tetrahedra_amplitude(tetrahedra, spin_data)
end
