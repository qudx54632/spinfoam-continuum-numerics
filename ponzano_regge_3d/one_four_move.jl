"""
Bookkeeping for the 1-4 Pachner move.

Before the move:

    (1, 2, 3, 4)

After inserting the internal vertex 0:

    (0, 1, 2, 3)
    (0, 1, 2, 4)
    (0, 1, 3, 4)
    (0, 2, 3, 4)

The boundary spin order follows the 6j convention for (1, 2, 3, 4):

    { j12  j13  j23 }
    { j34  j24  j14 }.
"""

const one_four_before_tetrahedron = (1, 2, 3, 4)

const one_four_after_tetrahedra = (
    (0, 1, 2, 3),
    (0, 1, 2, 4),
    (0, 1, 3, 4),
    (0, 2, 3, 4),
)

const one_four_boundary_edges = (
    (1, 2),
    (1, 3),
    (2, 3),
    (3, 4),
    (2, 4),
    (1, 4),
)

const one_four_internal_edges = (
    (0, 1),
    (0, 2),
    (0, 3),
    (0, 4),
)

"""Test admissibility of one 1-4 internal-spin assignment."""
function one_four_internal_spins_admissible(boundary_spins, internal_spins)
    j12, j13, j23, j34, j24, j14 = boundary_spins
    j01, j02, j03, j04 = internal_spins

    boundary_faces = (
        (j12, j13, j23),
        (j12, j14, j24),
        (j13, j14, j34),
        (j23, j24, j34),
    )

    internal_faces = (
        (j01, j02, j12),
        (j01, j03, j13),
        (j01, j04, j14),
        (j02, j03, j23),
        (j02, j04, j24),
        (j03, j04, j34),
    )

    return all(admissible_triangle, boundary_faces) &&
           all(admissible_triangle, internal_faces)
end

"""Build spin data from boundary spins and internal spins."""
function one_four_spin_data(boundary_spins, internal_spins)
    length(boundary_spins) == 6 || throw(ArgumentError("boundary_spins must have length 6"))
    length(internal_spins) == 4 || throw(ArgumentError("internal_spins must have length 4"))

    spin_data = Dict{Tuple{Int,Int},Real}()

    for (edge, spin) in zip(one_four_boundary_edges, boundary_spins)
        spin_data[edge] = spin
    end

    for (edge, spin) in zip(one_four_internal_edges, internal_spins)
        spin_data[edge] = spin
    end

    return spin_data
end

"""The before-move amplitude is the single boundary tetrahedron 6j symbol."""
function one_four_before_amplitude(boundary_spins)
    return tetrahedron_amplitude(boundary_spins[1:3], boundary_spins[4:6])
end

"""Fixed-spin after-move summand for the four tetrahedra around vertex 0."""
function one_four_after_amplitude(boundary_spins, internal_spins; chi0 = 0)
    spin_data = one_four_spin_data(boundary_spins, internal_spins)

    return state_sum_amplitude(
        one_four_after_tetrahedra,
        spin_data,
        one_four_internal_edges;
        chi0 = chi0,
    )
end

"""List admissible internal spins `(j01, j02, j03, j04)` up to `j_max`."""
function one_four_admissible_internal_spins(boundary_spins, j_max)
    spin_values = half_integer_spins(j_max)
    admissible_spins = NTuple{4,eltype(spin_values)}[]

    for j01 in spin_values, j02 in spin_values, j03 in spin_values, j04 in spin_values
        internal_spins = (j01, j02, j03, j04)

        if one_four_internal_spins_admissible(boundary_spins, internal_spins)
            push!(admissible_spins, internal_spins)
        end
    end

    return admissible_spins
end

"""Cutoff sum over admissible internal spins in the after-move triangulation."""
function one_four_after_sum(boundary_spins, j_max; chi0 = 0)
    total = 0

    for internal_spins in one_four_admissible_internal_spins(boundary_spins, j_max)
        total += one_four_after_amplitude(boundary_spins, internal_spins; chi0 = chi0)
    end

    return total
end

"""List admissible assignments with `j04` fixed and `j01,j02,j03 <= j_max`."""
function one_four_admissible_internal_spins_fixed_j04(boundary_spins, fixed_j04, j_max)
    spin_values = half_integer_spins(j_max)
    admissible_spins = []

    for j01 in spin_values, j02 in spin_values, j03 in spin_values
        internal_spins = (j01, j02, j03, fixed_j04)

        if one_four_internal_spins_admissible(boundary_spins, internal_spins)
            push!(admissible_spins, internal_spins)
        end
    end

    return admissible_spins
end

"""After-move sum with `j04` fixed."""
function one_four_after_sum_fixed_j04(boundary_spins, fixed_j04, j_max; chi0 = 0)
    total = 0

    for internal_spins in one_four_admissible_internal_spins_fixed_j04(boundary_spins, fixed_j04, j_max)
        total += one_four_after_amplitude(boundary_spins, internal_spins; chi0 = chi0)
    end

    return total
end

"""Gauge-volume-normalized amplitude from the fixed-`j04` sum."""
function one_four_physical_amplitude_fixed_j04(boundary_spins, fixed_j04, j_max; chi0 = 0)
    return one_four_after_sum_fixed_j04(boundary_spins, fixed_j04, j_max; chi0 = chi0) /
           edge_amplitude(fixed_j04)^2
end
