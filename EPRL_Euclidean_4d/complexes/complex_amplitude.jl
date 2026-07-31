"""Euclidean EPRL amplitude for a fixed labelled 4D complex."""

"""
Fixed-label amplitude of a multi-4-simplex complex.

The convention is

    A_C(j_f,i_tau)
      = prod_{f internal} d_{j_f}
        prod_{sigma in C} A_sigma^EPRL(j_f,i_tau).

Here `spin_data` assigns one diagonal SU(2) spin to every triangle appearing
in the complex, and `intertwiner_data` assigns one diagonal SU(2) intertwiner
to every tetrahedron appearing in the complex.
"""
function eprl_complex_amplitude(
    simplices,
    spin_data::AbstractDict,
    intertwiner_data::AbstractDict,
    internal_triangles,
    gamma,
)
    amplitude = triangle_amplitude_product(internal_triangles, spin_data)

    for simplex in simplices
        vertex = eprl_vertex_amplitude(simplex, spin_data, intertwiner_data, gamma)
        vertex == 0.0 && return 0.0

        amplitude *= vertex
    end

    return amplitude
end

"""
Sum the fixed-spin complex amplitude over internal diagonal intertwiners.

This is useful after the internal triangle spins have already been chosen.
The remaining state-sum variables are the intertwiners on internal tetrahedra.
"""
function eprl_complex_sum_over_internal_intertwiners(
    simplices,
    spin_data::AbstractDict,
    boundary_intertwiner_data::AbstractDict,
    internal_tetrahedra,
    internal_triangles,
    gamma,
)
    intertwiner_ranges = [
        allowed_global_intertwiners(tetrahedron, spin_data)
        for tetrahedron in internal_tetrahedra
    ]

    any(isempty, intertwiner_ranges) && return 0.0

    total = 0.0

    for internal_intertwiners in Iterators.product(intertwiner_ranges...)
        intertwiner_data = copy(boundary_intertwiner_data)

        for (tetrahedron, i) in zip(internal_tetrahedra, internal_intertwiners)
            intertwiner_data[tetrahedron] = i
        end

        total += eprl_complex_amplitude(
            simplices,
            spin_data,
            intertwiner_data,
            internal_triangles,
            gamma,
        )
    end

    return total
end
