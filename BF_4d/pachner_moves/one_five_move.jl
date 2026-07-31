"""
1-5 Pachner move for the 4D SU(2) BF state sum.

The after side has ten internal triangle spins.  In the gauge-fixed check
used here, the four redundant triangles containing the internal edge (0,1)
are fixed, and the remaining six internal triangle spins are summed up to a
cutoff.
"""

# Oriented 4-simplices, using the convention in the note.
const one_five_before_simplex_oriented = (1, 2, 3, 4, 5)

const one_five_after_simplices_oriented = (
    (0, 2, 3, 4, 5),
    (1, 0, 3, 4, 5),
    (0, 1, 2, 4, 5),
    (1, 0, 2, 3, 5),
    (0, 1, 2, 3, 4),
)

# Boundary spin input order: local spin order of the before 4-simplex.
const one_five_boundary_triangles = (
    (3, 4, 5),
    (2, 4, 5),
    (2, 3, 5),
    (2, 3, 4),
    (1, 4, 5),
    (1, 3, 5),
    (1, 3, 4),
    (1, 2, 5),
    (1, 2, 4),
    (1, 2, 3),
)

# Boundary intertwiner input order: tetrahedra opposite vertices 1,...,5.
const one_five_boundary_tetrahedra = (
    (2, 3, 4, 5),
    (1, 3, 4, 5),
    (1, 2, 4, 5),
    (1, 2, 3, 5),
    (1, 2, 3, 4),
)

# Gauge-fixed test: these four redundant internal triangle spins are fixed.
const one_five_fixed_internal_triangles = (
    (0, 1, 2),
    (0, 1, 3),
    (0, 1, 4),
    (0, 1, 5),
)

# The remaining six internal triangle spins are summed.
const one_five_free_internal_triangles = (
    (0, 2, 3),
    (0, 2, 4),
    (0, 2, 5),
    (0, 3, 4),
    (0, 3, 5),
    (0, 4, 5),
)

const one_five_internal_triangles = (
    one_five_fixed_internal_triangles...,
    one_five_free_internal_triangles...,
)

const one_five_internal_tetrahedra = (
    (0, 1, 2, 3),
    (0, 1, 2, 4),
    (0, 1, 2, 5),
    (0, 1, 3, 4),
    (0, 1, 3, 5),
    (0, 1, 4, 5),
    (0, 2, 3, 4),
    (0, 2, 3, 5),
    (0, 2, 4, 5),
    (0, 3, 4, 5),
)

function one_five_boundary_spin_data(boundary_spins)
    length(boundary_spins) == length(one_five_boundary_triangles) ||
        throw(ArgumentError("boundary_spins must have length 10"))

    return Dict{Tuple{Int, Int, Int}, Real}(
        triangle => spin
        for (triangle, spin) in zip(one_five_boundary_triangles, boundary_spins)
    )
end

function one_five_boundary_intertwiner_data(boundary_intertwiners)
    length(boundary_intertwiners) == length(one_five_boundary_tetrahedra) ||
        throw(ArgumentError("boundary_intertwiners must have length 5"))

    return Dict{Tuple{Int, Int, Int, Int}, Real}(
        tetrahedron => i
        for (tetrahedron, i) in zip(one_five_boundary_tetrahedra, boundary_intertwiners)
    )
end

"""Insert the ten after-move internal triangle spins into the boundary data."""
function one_five_spin_data_after(boundary_spins, internal_spins)
    length(internal_spins) == length(one_five_internal_triangles) ||
        throw(ArgumentError("internal_spins must have length 10"))

    spin_data = one_five_boundary_spin_data(boundary_spins)

    for (triangle, spin) in zip(one_five_internal_triangles, internal_spins)
        spin_data[triangle] = spin
    end

    return spin_data
end

"""Insert the ten after-move internal intertwiners into the boundary data."""
function one_five_intertwiner_data_after(boundary_intertwiners, internal_intertwiners)
    length(internal_intertwiners) == length(one_five_internal_tetrahedra) ||
        throw(ArgumentError("internal_intertwiners must have length 10"))

    intertwiner_data = one_five_boundary_intertwiner_data(boundary_intertwiners)

    for (tetrahedron, i) in zip(one_five_internal_tetrahedra, internal_intertwiners)
        intertwiner_data[tetrahedron] = i
    end

    return intertwiner_data
end

"""Before amplitude: one oriented 4-simplex and no internal labels."""
function one_five_before_amplitude_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners;
    state_sum_function = state_sum_amplitude_recoupled_oriented,
)
    spin_data = one_five_boundary_spin_data(boundary_spins)
    intertwiner_data = one_five_boundary_intertwiner_data(boundary_intertwiners)

    return state_sum_function(
        (one_five_before_simplex_oriented,),
        spin_data,
        intertwiner_data,
        (),
    )
end

"""Put the four fixed and six free internal triangle spins into one tuple."""
function one_five_internal_spins(fixed_internal_spins, free_internal_spins)
    length(fixed_internal_spins) == 4 ||
        throw(ArgumentError("fixed_internal_spins must have length 4"))
    length(free_internal_spins) == 6 ||
        throw(ArgumentError("free_internal_spins must have length 6"))

    return (fixed_internal_spins..., free_internal_spins...)
end

"""
Gauge-fixed after amplitude.

The first four internal triangle spins in `one_five_internal_triangles` are
fixed.  The remaining six are summed over half-integers up to `j_max`; for
each admissible spin assignment, the ten internal tetrahedron intertwiners are
summed.
"""
function one_five_after_sum_fixed_internal_triangles_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    fixed_internal_spins,
    j_max;
    state_sum_function = state_sum_amplitude_recoupled_oriented,
)
    spin_values = half_integer_spins(j_max)
    free_spin_ranges = ntuple(_ -> spin_values, 6)
    total = 0.0

    for free_internal_spins in Iterators.product(free_spin_ranges...)
        internal_spins = one_five_internal_spins(
            fixed_internal_spins,
            free_internal_spins,
        )

        spin_data = one_five_spin_data_after(boundary_spins, internal_spins)

        i_ranges = [
            allowed_global_intertwiners(tetrahedron, spin_data)
            for tetrahedron in one_five_internal_tetrahedra
        ]

        any(isempty, i_ranges) && continue

        for internal_intertwiners in Iterators.product(i_ranges...)
            intertwiner_data = one_five_intertwiner_data_after(
                boundary_intertwiners,
                internal_intertwiners,
            )

            total += state_sum_function(
                one_five_after_simplices_oriented,
                spin_data,
                intertwiner_data,
                one_five_internal_triangles,
            )
        end
    end

    return total
end

"""
Fixed-spin factor for the four redundant internal triangle spins.

Each fixed Fourier component contributes `d_J chi_J(I) = d_J^2`.
"""
function one_five_fixed_internal_triangle_factor(fixed_internal_spins)
    return prod(dimension(j)^2 for j in fixed_internal_spins)
end
