"""1-5 Pachner move for the Euclidean EPRL state sum."""

# The before side is the same single oriented 4-simplex used in main.jl.
const one_five_before_simplex = (1, 2, 3, 4, 5)

# Oriented cone subdivision of the same boundary 4-simplex.
const one_five_after_simplices = (
    (0, 2, 3, 4, 5),
    (1, 0, 3, 4, 5),
    (0, 1, 2, 4, 5),
    (1, 0, 2, 3, 5),
    (0, 1, 2, 3, 4),
)

# Boundary labels are given in the local order of the before 4-simplex.
const one_five_boundary_triangles = oriented_simplex_spin_keys(one_five_before_simplex)
const one_five_boundary_tetrahedra = oriented_simplex_tetrahedra(one_five_before_simplex)

# Gauge-fixed comparison: fix the four internal triangles containing the
# internal edge (0,1).
const one_five_fixed_internal_triangles = (
    (0, 1, 2),
    (0, 1, 3),
    (0, 1, 4),
    (0, 1, 5),
)

# The remaining six internal triangle spins are summed up to a cutoff.
const one_five_free_internal_triangles = (
    (0, 2, 3),
    (0, 2, 4),
    (0, 2, 5),
    (0, 3, 4),
    (0, 3, 5),
    (0, 4, 5),
)

# All internal triangles of the after move.
const one_five_internal_triangles = (
    one_five_fixed_internal_triangles...,
    one_five_free_internal_triangles...,
)

# Internal tetrahedra are the ten tetrahedra containing 0.
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

    return Dict(
        triangle => spin
        for (triangle, spin) in zip(one_five_boundary_triangles, boundary_spins)
    )
end

function one_five_boundary_intertwiner_data(boundary_intertwiners)
    length(boundary_intertwiners) == length(one_five_boundary_tetrahedra) ||
        throw(ArgumentError("boundary_intertwiners must have length 5"))

    return Dict(
        tetrahedron => i
        for (tetrahedron, i) in zip(one_five_boundary_tetrahedra, boundary_intertwiners)
    )
end

function one_five_after_spin_data(boundary_spins, internal_spins)
    length(internal_spins) == length(one_five_internal_triangles) ||
        throw(ArgumentError("internal_spins must have length 10"))

    spin_data = one_five_boundary_spin_data(boundary_spins)

    for (triangle, spin) in zip(one_five_internal_triangles, internal_spins)
        spin_data[triangle] = spin
    end

    return spin_data
end

function one_five_internal_spins(fixed_internal_spins, free_internal_spins)
    length(fixed_internal_spins) == length(one_five_fixed_internal_triangles) ||
        throw(ArgumentError("fixed_internal_spins must have length 4"))
    length(free_internal_spins) == length(one_five_free_internal_triangles) ||
        throw(ArgumentError("free_internal_spins must have length 6"))

    return (fixed_internal_spins..., free_internal_spins...)
end

"""Allowed diagonal internal spins below the cutoff for this gamma."""
function one_five_internal_spin_values(j_max, gamma)
    return [
        j
        for j in half_integer_spins(j_max)
        if eprl_spin_allowed(j, gamma)
    ]
end

"""Before amplitude: the original single 4-simplex."""
function one_five_before_amplitude(boundary_spins, boundary_intertwiners, gamma)
    spin_data = one_five_boundary_spin_data(boundary_spins)
    intertwiner_data = one_five_boundary_intertwiner_data(boundary_intertwiners)

    return eprl_complex_amplitude(
        (one_five_before_simplex,),
        spin_data,
        intertwiner_data,
        (),
        gamma,
    )
end

"""
After amplitude with an internal spin cutoff.

The boundary data are the same as on the before 4-simplex.  The ten internal
triangle spins are summed up to `j_max`, and the ten internal tetrahedron
intertwiners are summed exactly for each spin assignment.
"""
function one_five_after_amplitude(boundary_spins, boundary_intertwiners, j_max, gamma)
    spin_values = one_five_internal_spin_values(j_max, gamma)
    spin_ranges = ntuple(_ -> spin_values, length(one_five_internal_triangles))

    boundary_intertwiner_data =
        one_five_boundary_intertwiner_data(boundary_intertwiners)

    total = 0.0

    for internal_spins in Iterators.product(spin_ranges...)
        spin_data = one_five_after_spin_data(boundary_spins, internal_spins)

        total += eprl_complex_sum_over_internal_intertwiners(
            one_five_after_simplices,
            spin_data,
            boundary_intertwiner_data,
            one_five_internal_tetrahedra,
            one_five_internal_triangles,
            gamma,
        )
    end

    return total
end

"""
After amplitude with four fixed internal triangle spins.

The fixed spins are assigned to

    (012), (013), (014), (015).

The six remaining internal triangle spins are summed up to `j_max`, and the
ten internal tetrahedron intertwiners are summed exactly for each spin
assignment.
"""
function one_five_after_amplitude_fixed_internal_spins(
    boundary_spins,
    boundary_intertwiners,
    fixed_internal_spins,
    j_max,
    gamma,
)
    all(eprl_spin_allowed(j, gamma) for j in fixed_internal_spins) ||
        return 0.0

    spin_values = one_five_internal_spin_values(j_max, gamma)
    free_spin_ranges = ntuple(_ -> spin_values, length(one_five_free_internal_triangles))

    boundary_intertwiner_data =
        one_five_boundary_intertwiner_data(boundary_intertwiners)

    total = 0.0

    for free_internal_spins in Iterators.product(free_spin_ranges...)
        internal_spins = one_five_internal_spins(
            fixed_internal_spins,
            free_internal_spins,
        )

        spin_data = one_five_after_spin_data(boundary_spins, internal_spins)

        total += eprl_complex_sum_over_internal_intertwiners(
            one_five_after_simplices,
            spin_data,
            boundary_intertwiner_data,
            one_five_internal_tetrahedra,
            one_five_internal_triangles,
            gamma,
        )
    end

    return total
end

"""
Cutoff scan for the fixed-internal-spin after amplitude.

Each free internal spin assignment is evaluated only once.  Its contribution
enters every cutoff larger than or equal to the largest free spin in that
assignment.
"""
function one_five_after_fixed_internal_spin_cutoff_scan(
    boundary_spins,
    boundary_intertwiners,
    fixed_internal_spins,
    cutoffs,
    gamma,
)
    all(eprl_spin_allowed(j, gamma) for j in fixed_internal_spins) ||
        return [(cutoff = cutoff, amplitude = 0.0) for cutoff in cutoffs]

    cutoffs = unique(sort(collect(cutoffs)))
    max_cutoff = maximum(cutoffs)
    spin_values = one_five_internal_spin_values(max_cutoff, gamma)
    free_spin_ranges = ntuple(_ -> spin_values, length(one_five_free_internal_triangles))

    boundary_intertwiner_data =
        one_five_boundary_intertwiner_data(boundary_intertwiners)

    shell_contribution = Dict(spin => 0.0 for spin in spin_values)

    for free_internal_spins in Iterators.product(free_spin_ranges...)
        internal_spins = one_five_internal_spins(
            fixed_internal_spins,
            free_internal_spins,
        )

        spin_data = one_five_after_spin_data(boundary_spins, internal_spins)

        contribution = eprl_complex_sum_over_internal_intertwiners(
            one_five_after_simplices,
            spin_data,
            boundary_intertwiner_data,
            one_five_internal_tetrahedra,
            one_five_internal_triangles,
            gamma,
        )

        shell = maximum(free_internal_spins)
        shell_contribution[shell] += contribution
    end

    return [
        (
            cutoff = cutoff,
            amplitude = sum(
                contribution
                for (shell, contribution) in shell_contribution
                if shell <= cutoff;
                init = 0.0,
            ),
        )
        for cutoff in cutoffs
    ]
end
