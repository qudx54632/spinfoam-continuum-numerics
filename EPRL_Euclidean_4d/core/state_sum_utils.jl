"""Return a triangle in canonical vertex order."""
function triangle_key(a::Integer, b::Integer, c::Integer)
    length(unique((a, b, c))) == 3 ||
        throw(ArgumentError("a triangle needs three different vertices"))

    return Tuple(sort([a, b, c]))
end

"""Return a tetrahedron in canonical vertex order."""
function tetrahedron_key(a::Integer, b::Integer, c::Integer, d::Integer)
    length(unique((a, b, c, d))) == 4 ||
        throw(ArgumentError("a tetrahedron needs four different vertices"))

    return Tuple(sort([a, b, c, d]))
end

"""Return a 4-simplex in canonical vertex order."""
function four_simplex_key(simplex)
    length(simplex) == 5 || throw(ArgumentError("a 4-simplex needs five vertices"))
    length(unique(simplex)) == 5 ||
        throw(ArgumentError("a 4-simplex needs five different vertices"))

    return Tuple(sort(collect(simplex)))
end

"""Multiply triangle dimensions `d_j = 2j + 1` over selected triangles."""
function triangle_amplitude_product(triangles, spin_data::AbstractDict)
    return prod((dimension(spin_data[triangle]) for triangle in triangles); init = 1)
end
