"""
2-4 Pachner move for the 4D SU(2) BF state sum.

The common files have already defined:

    state_sum_amplitude_recoupled_oriented(...)
    allowed_global_intertwiners(...)

This file only specifies the 2-4 geometry and performs the two sums.
"""

# Oriented 4-simplices, using the convention in the note.
const two_four_before_simplices_oriented = (
    (1, 0, 2, 3, 4),
    (0, 1, 2, 3, 5),
)

const two_four_after_simplices_oriented = (
    (2, 1, 3, 4, 5),
    (0, 2, 3, 4, 5),
    (1, 0, 3, 4, 5),
    (0, 1, 2, 4, 5),
)

# Boundary spin input order.
const two_four_boundary_triangles = (
    (0, 1, 2),
    (0, 1, 3),
    (0, 2, 3),
    (1, 2, 3),
    (0, 1, 4),
    (0, 2, 4),
    (0, 3, 4),
    (1, 2, 4),
    (1, 3, 4),
    (2, 3, 4),
    (0, 1, 5),
    (0, 2, 5),
    (0, 3, 5),
    (1, 2, 5),
    (1, 3, 5),
    (2, 3, 5),
)

# Boundary intertwiner input order.
const two_four_boundary_tetrahedra = (
    (0, 1, 2, 4),
    (0, 1, 3, 4),
    (0, 2, 3, 4),
    (1, 2, 3, 4),
    (0, 1, 2, 5),
    (0, 1, 3, 5),
    (0, 2, 3, 5),
    (1, 2, 3, 5),
)

# Before side: one internal tetrahedron and no internal triangle.
const two_four_before_internal_triangles = ()
const two_four_before_internal_tetrahedron = (0, 1, 2, 3)

# After side: four internal triangles around the internal edge (4,5),
# and six internal tetrahedra.
const two_four_after_internal_triangles = (
    (0, 4, 5),
    (1, 4, 5),
    (2, 4, 5),
    (3, 4, 5),
)

const two_four_after_internal_tetrahedra = (
    (0, 1, 4, 5),
    (0, 2, 4, 5),
    (0, 3, 4, 5),
    (1, 2, 4, 5),
    (1, 3, 4, 5),
    (2, 3, 4, 5),
)

function two_four_boundary_spin_data(boundary_spins)
    length(boundary_spins) == length(two_four_boundary_triangles) ||
        throw(ArgumentError("boundary_spins must have length 16"))

    return Dict{Tuple{Int, Int, Int}, Real}(
        triangle => spin
        for (triangle, spin) in zip(two_four_boundary_triangles, boundary_spins)
    )
end

function two_four_boundary_intertwiner_data(boundary_intertwiners)
    length(boundary_intertwiners) == length(two_four_boundary_tetrahedra) ||
        throw(ArgumentError("boundary_intertwiners must have length 8"))

    return Dict{Tuple{Int, Int, Int, Int}, Real}(
        tetrahedron => i
        for (tetrahedron, i) in zip(two_four_boundary_tetrahedra, boundary_intertwiners)
    )
end

"""Insert the after-move internal triangle spins into the boundary spin data."""
function two_four_spin_data_after(boundary_spins, internal_spins)
    length(internal_spins) == length(two_four_after_internal_triangles) ||
        throw(ArgumentError("internal_spins must have length 4"))

    spin_data = two_four_boundary_spin_data(boundary_spins)

    for (triangle, spin) in zip(two_four_after_internal_triangles, internal_spins)
        spin_data[triangle] = spin
    end

    return spin_data
end

"""Insert the after-move internal intertwiners into the boundary data."""
function two_four_intertwiner_data_after(boundary_intertwiners, internal_intertwiners)
    length(internal_intertwiners) == length(two_four_after_internal_tetrahedra) ||
        throw(ArgumentError("internal_intertwiners must have length 6"))

    intertwiner_data = two_four_boundary_intertwiner_data(boundary_intertwiners)

    for (tetrahedron, i) in zip(two_four_after_internal_tetrahedra, internal_intertwiners)
        intertwiner_data[tetrahedron] = i
    end

    return intertwiner_data
end

"""
Before amplitude.

Only the internal tetrahedron intertwiner i_0123 is summed.
"""
function two_four_before_sum_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners;
    state_sum_function = state_sum_amplitude_recoupled_oriented,
)
    spin_data = two_four_boundary_spin_data(boundary_spins)
    i_values = allowed_global_intertwiners(two_four_before_internal_tetrahedron, spin_data)
    total = 0.0

    for i0123 in i_values
        intertwiner_data = two_four_boundary_intertwiner_data(boundary_intertwiners)
        intertwiner_data[two_four_before_internal_tetrahedron] = i0123

        total += state_sum_function(
            two_four_before_simplices_oriented,
            spin_data,
            intertwiner_data,
            two_four_before_internal_triangles,
        )
    end

    return total
end

"""Put `fixed_spin` into one of the four after-move internal triangle slots."""
function two_four_internal_spins_with_fixed_triangle(fixed_position, fixed_spin, free_spins)
    fixed_position in 1:4 ||
        throw(ArgumentError("fixed_position must be between 1 and 4"))
    length(free_spins) == 3 ||
        throw(ArgumentError("free_spins must have length 3"))

    internal_spins = Vector{Real}(undef, 4)
    free_position = 1

    for position in 1:4
        if position == fixed_position
            internal_spins[position] = fixed_spin
        else
            internal_spins[position] = free_spins[free_position]
            free_position += 1
        end
    end

    return Tuple(internal_spins)
end

"""
Gauge-fixed after amplitude.

One internal triangle spin is fixed.  The other three internal triangle spins
are summed over half-integers up to `j_max`; for each admissible spin
assignment, the six internal tetrahedron intertwiners are summed.
"""
function two_four_after_sum_fixed_internal_triangle_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    fixed_position,
    fixed_spin,
    j_max;
    state_sum_function = state_sum_amplitude_recoupled_oriented,
)
    spin_values = half_integer_spins(j_max)
    total = 0.0

    for free_spins in Iterators.product(spin_values, spin_values, spin_values)
        internal_spins = two_four_internal_spins_with_fixed_triangle(
            fixed_position,
            fixed_spin,
            free_spins,
        )

        spin_data = two_four_spin_data_after(boundary_spins, internal_spins)

        i_ranges = [
            allowed_global_intertwiners(tetrahedron, spin_data)
            for tetrahedron in two_four_after_internal_tetrahedra
        ]

        any(isempty, i_ranges) && continue

        for internal_intertwiners in Iterators.product(i_ranges...)
            intertwiner_data = two_four_intertwiner_data_after(
                boundary_intertwiners,
                internal_intertwiners,
            )

            total += state_sum_function(
                two_four_after_simplices_oriented,
                spin_data,
                intertwiner_data,
                two_four_after_internal_triangles,
            )
        end
    end

    return total
end
