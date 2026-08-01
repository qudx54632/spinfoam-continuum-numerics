"""
Livine-Speziale coherent boundary states in the fixed global basis.

Reference:
    E. R. Livine and S. Speziale, arXiv:0705.0674.

The code uses the highest-weight SU(2) coherent state

    |j,n> = g(n) |j,j>,

where `g(n)` rotates the z-axis to the unit vector `n`.  In the magnetic
basis this gives

    <j,m|j,n>
      = sqrt((2j)!/((j+m)!(j-m)!))
        cos(theta/2)^(j+m) sin(theta/2)^(j-m)
        exp[-i (j-m) phi].

The coherent intertwiner coefficient used in the vertex amplitude is

    c_i(j_l,n_l)
      = sum_{m_l} conj(i_{m_1 m_2 m_3 m_4})
                  prod_l <j_l,m_l|j_l,n_l>,

with the magnetic constraint `m1+m2+m3+m4=0`.  The intertwiner channel is
the normalized four-valent channel `(j1,j2) -> i -> (j3,j4)`.
"""

const coherent_state_coefficient_cache = Dict()
const coherent_intertwiner_coefficient_cache = Dict()

empty!(coherent_state_coefficient_cache)
empty!(coherent_intertwiner_coefficient_cache)

function spin_integer(x, name)
    isinteger(x) || throw(ArgumentError("$name must be an integer"))
    return Int(x)
end

function unit_vector3(n)
    length(n) == 3 || throw(ArgumentError("a normal vector must have length 3"))

    x, y, z = Float64.(n)
    norm = sqrt(x^2 + y^2 + z^2)
    norm > 0 || throw(ArgumentError("a normal vector must be non-zero"))

    return (x / norm, y / norm, z / norm)
end

function spherical_angles(n)
    x, y, z = unit_vector3(n)
    theta = acos(clamp(z, -1.0, 1.0))
    phi = atan(y, x)

    return theta, phi
end

"""Coefficient `<j,m|j,n>` of the SU(2) coherent state."""
function coherent_state_coefficient(j, m, n)
    allowed_magnetic_label(j, m) ||
        throw(ArgumentError("m is not an allowed magnetic label for j"))

    normal = unit_vector3(n)
    key = (j, m, normal)

    if haskey(coherent_state_coefficient_cache, key)
        return coherent_state_coefficient_cache[key]
    end

    theta, phi = spherical_angles(normal)
    two_j = spin_integer(2j, "2j")
    j_plus_m = spin_integer(j + m, "j+m")
    j_minus_m = spin_integer(j - m, "j-m")

    value =
        sqrt(Float64(binomial(two_j, j_plus_m))) *
        cos(theta / 2)^j_plus_m *
        sin(theta / 2)^j_minus_m *
        exp(-im * j_minus_m * phi)

    coherent_state_coefficient_cache[key] = value

    return value
end

"""Coefficient `c_i(j_l,n_l)` of the coherent intertwiner in channel `i`."""
function coherent_intertwiner_coefficient(spins, i, normals)
    spins = Tuple(spins)
    normals = Tuple(normals)

    length(spins) == 4 || throw(ArgumentError("spins must have length 4"))
    length(normals) == 4 || throw(ArgumentError("normals must have length 4"))

    admissible_4valent(spins, i) || return 0.0 + 0.0im

    normalized_normals = Tuple(unit_vector3(n) for n in normals)
    key = (spins, i, normalized_normals)

    if haskey(coherent_intertwiner_coefficient_cache, key)
        return coherent_intertwiner_coefficient_cache[key]
    end

    j1, j2, j3, j4 = spins
    total = 0.0 + 0.0im

    for m1 in magnetic_values(j1),
        m2 in magnetic_values(j2),
        m3 in magnetic_values(j3)

        m4 = -m1 - m2 - m3
        allowed_magnetic_label(j4, m4) || continue

        magnetic_labels = (m1, m2, m3, m4)
        intertwiner = four_valent_intertwiner(spins, i, magnetic_labels)
        intertwiner == 0.0 && continue

        coherent_product = prod(
            coherent_state_coefficient(j, m, n)
            for (j, m, n) in zip(spins, magnetic_labels, normalized_normals);
            init = 1.0 + 0.0im,
        )

        total += conj(intertwiner) * coherent_product
    end

    coherent_intertwiner_coefficient_cache[key] = total

    return total
end

function tetrahedron_triangle_normal(normal_data::AbstractDict, tetrahedron, triangle)
    tau = tetrahedron_key(tetrahedron...)
    f = triangle_key(triangle...)

    if haskey(normal_data, (tau, f))
        return normal_data[(tau, f)]
    elseif haskey(normal_data, f)
        return normal_data[f]
    else
        throw(KeyError((tau, f)))
    end
end

"""
Same coefficient, with spins and normals read in the fixed global tetrahedron
basis.

The preferred normal-data key is `(tetrahedron, triangle)`, because the same
boundary triangle has different normals when seen from its two neighbouring
tetrahedra.
"""
function coherent_intertwiner_coefficient(
    tetrahedron,
    spin_data::AbstractDict,
    normal_data::AbstractDict,
    i,
)
    spin_keys = global_intertwiner_spin_keys(tetrahedron)
    spins = Tuple(spin_data[key] for key in spin_keys)
    normals = Tuple(
        tetrahedron_triangle_normal(normal_data, tetrahedron, key)
        for key in spin_keys
    )

    return coherent_intertwiner_coefficient(spins, i, normals)
end

"""All coherent-intertwiner coefficients allowed by the tetrahedron spins."""
function coherent_intertwiner_coefficients(tetrahedron, spin_data::AbstractDict, normal_data::AbstractDict)
    return [
        (
            i = i,
            coefficient = coherent_intertwiner_coefficient(
                tetrahedron,
                spin_data,
                normal_data,
                i,
            ),
        )
        for i in allowed_global_intertwiners(tetrahedron, spin_data)
    ]
end

"""Closure vector `sum_l j_l n_l` for one coherent tetrahedron."""
function closure_vector(spins, normals)
    spins = Tuple(spins)
    normals = Tuple(unit_vector3(n) for n in normals)

    length(spins) == 4 || throw(ArgumentError("spins must have length 4"))
    length(normals) == 4 || throw(ArgumentError("normals must have length 4"))

    x = sum(Float64(j) * n[1] for (j, n) in zip(spins, normals); init = 0.0)
    y = sum(Float64(j) * n[2] for (j, n) in zip(spins, normals); init = 0.0)
    z = sum(Float64(j) * n[3] for (j, n) in zip(spins, normals); init = 0.0)

    return (x, y, z)
end
