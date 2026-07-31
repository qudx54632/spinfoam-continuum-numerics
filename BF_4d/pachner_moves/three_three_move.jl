"""
3-3 Pachner move for the 4D SU(2) BF state sum.

Both sides have one internal triangle spin and three internal tetrahedron
intertwiners.  This file only writes the 3-3 geometry and performs those
finite sums.
"""

# Oriented 4-simplices, using the convention in the note.
const three_three_before_simplices_oriented = (
    (1, 0, 2, 3, 4),
    (0, 1, 2, 3, 5),
    (1, 0, 2, 4, 5),
)

const three_three_after_simplices_oriented = (
    (2, 1, 3, 4, 5),
    (0, 2, 3, 4, 5),
    (1, 0, 3, 4, 5),
)

# Boundary spin input order.
const three_three_boundary_triangles = (
    (0, 1, 3),
    (0, 1, 4),
    (0, 1, 5),
    (0, 2, 3),
    (0, 2, 4),
    (0, 2, 5),
    (1, 2, 3),
    (1, 2, 4),
    (1, 2, 5),
    (0, 3, 4),
    (0, 3, 5),
    (0, 4, 5),
    (1, 3, 4),
    (1, 3, 5),
    (1, 4, 5),
    (2, 3, 4),
    (2, 3, 5),
    (2, 4, 5),
)

# Boundary intertwiner input order.
const three_three_boundary_tetrahedra = (
    (0, 1, 3, 4),
    (0, 1, 3, 5),
    (0, 1, 4, 5),
    (0, 2, 3, 4),
    (0, 2, 3, 5),
    (0, 2, 4, 5),
    (1, 2, 3, 4),
    (1, 2, 3, 5),
    (1, 2, 4, 5),
)

# Before side: internal triangle (0,1,2).
const three_three_before_internal_triangles = ((0, 1, 2),)
const three_three_before_internal_tetrahedra = (
    (0, 1, 2, 3),
    (0, 1, 2, 4),
    (0, 1, 2, 5),
)

# After side: internal triangle (3,4,5).
const three_three_after_internal_triangles = ((3, 4, 5),)
const three_three_after_internal_tetrahedra = (
    (0, 3, 4, 5),
    (1, 3, 4, 5),
    (2, 3, 4, 5),
)

function three_three_boundary_spin_data(boundary_spins)
    length(boundary_spins) == length(three_three_boundary_triangles) ||
        throw(ArgumentError("boundary_spins must have length 18"))

    return Dict{Tuple{Int, Int, Int}, Real}(
        triangle => spin
        for (triangle, spin) in zip(three_three_boundary_triangles, boundary_spins)
    )
end

function three_three_boundary_intertwiner_data(boundary_intertwiners)
    length(boundary_intertwiners) == length(three_three_boundary_tetrahedra) ||
        throw(ArgumentError("boundary_intertwiners must have length 9"))

    return Dict{Tuple{Int, Int, Int, Int}, Real}(
        tetrahedron => i
        for (tetrahedron, i) in zip(three_three_boundary_tetrahedra, boundary_intertwiners)
    )
end

"""Insert one internal triangle spin into the boundary spin data."""
function three_three_spin_data(boundary_spins, internal_triangle, internal_spin)
    spin_data = three_three_boundary_spin_data(boundary_spins)
    spin_data[internal_triangle] = internal_spin

    return spin_data
end

"""Insert the three internal tetrahedron intertwiners into the boundary data."""
function three_three_intertwiner_data(
    boundary_intertwiners,
    internal_tetrahedra,
    internal_intertwiners,
)
    length(internal_intertwiners) == length(internal_tetrahedra) ||
        throw(ArgumentError("wrong number of internal intertwiners"))

    intertwiner_data = three_three_boundary_intertwiner_data(boundary_intertwiners)

    for (tetrahedron, i) in zip(internal_tetrahedra, internal_intertwiners)
        intertwiner_data[tetrahedron] = i
    end

    return intertwiner_data
end

"""
Before amplitude.

Sum over the internal triangle spin j_012 and over the three internal
tetrahedron intertwiners.
"""
function three_three_before_sum_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    j_max;
    state_sum_function = state_sum_amplitude_recoupled_oriented,
)
    total = 0.0

    for j012 in half_integer_spins(j_max)
        spin_data = three_three_spin_data(boundary_spins, (0, 1, 2), j012)

        i_ranges = [
            allowed_global_intertwiners(tetrahedron, spin_data)
            for tetrahedron in three_three_before_internal_tetrahedra
        ]

        any(isempty, i_ranges) && continue

        for internal_intertwiners in Iterators.product(i_ranges...)
            intertwiner_data = three_three_intertwiner_data(
                boundary_intertwiners,
                three_three_before_internal_tetrahedra,
                internal_intertwiners,
            )

            total += state_sum_function(
                three_three_before_simplices_oriented,
                spin_data,
                intertwiner_data,
                three_three_before_internal_triangles,
            )
        end
    end

    return total
end

"""
After amplitude.

Sum over the internal triangle spin j_345 and over the three internal
tetrahedron intertwiners.
"""
function three_three_after_sum_oriented_recoupled(
    boundary_spins,
    boundary_intertwiners,
    j_max;
    state_sum_function = state_sum_amplitude_recoupled_oriented,
)
    total = 0.0

    for j345 in half_integer_spins(j_max)
        spin_data = three_three_spin_data(boundary_spins, (3, 4, 5), j345)

        i_ranges = [
            allowed_global_intertwiners(tetrahedron, spin_data)
            for tetrahedron in three_three_after_internal_tetrahedra
        ]

        any(isempty, i_ranges) && continue

        for internal_intertwiners in Iterators.product(i_ranges...)
            intertwiner_data = three_three_intertwiner_data(
                boundary_intertwiners,
                three_three_after_internal_tetrahedra,
                internal_intertwiners,
            )

            total += state_sum_function(
                three_three_after_simplices_oriented,
                spin_data,
                intertwiner_data,
                three_three_after_internal_triangles,
            )
        end
    end

    return total
end
