const simplex_amplitude_recoupled_oriented_6j_cache = Dict()

empty!(simplex_amplitude_recoupled_oriented_6j_cache)

include(joinpath(@__DIR__, "recoupling_6j.jl"))

"""
Orientation-aware BF 4-simplex vertex in the fixed global tetrahedron basis,
using the fast local 6j expression.

This has the same input/output convention as
`simplex_amplitude_recoupled_oriented`; only the local 15j evaluation is
changed.
"""
function simplex_amplitude_recoupled_oriented_6j(
    simplex,
    spin_data::AbstractDict,
    global_intertwiner_data::AbstractDict,
)
    sigma = oriented_four_simplex(simplex)
    tetrahedra = oriented_simplex_tetrahedra(sigma)

    J = oriented_simplex_spins(sigma, spin_data)
    i_tau = Tuple(global_intertwiner_data[tau] for tau in tetrahedra)

    cache_key = (sigma, J, i_tau)
    if haskey(simplex_amplitude_recoupled_oriented_6j_cache, cache_key)
        return simplex_amplitude_recoupled_oriented_6j_cache[cache_key]
    end

    k_ranges = [
        allowed_oriented_local_intertwiners(sigma, tau, spin_data)
        for tau in tetrahedra
    ]

    if any(isempty, k_ranges)
        simplex_amplitude_recoupled_oriented_6j_cache[cache_key] = 0.0
        return 0.0
    end

    amplitude = 0.0

    for k_tau in Iterators.product(k_ranges...)
        R = recoupling_product_6j(sigma, tetrahedra, spin_data, i_tau, k_tau)
        R == 0.0 && continue

        amplitude += R * local_vertex_amplitude_6j(J, k_tau)
    end

    simplex_amplitude_recoupled_oriented_6j_cache[cache_key] = amplitude

    return amplitude
end

"""Fixed-label BF state-sum summand using recoupled oriented 6j vertices."""
function state_sum_amplitude_recoupled_oriented_6j(
    simplices,
    spin_data::AbstractDict,
    global_intertwiner_data::AbstractDict,
    internal_triangles,
)
    amplitude = triangle_amplitude_product(internal_triangles, spin_data)

    for simplex in simplices
        amplitude *= simplex_amplitude_recoupled_oriented_6j(
            simplex,
            spin_data,
            global_intertwiner_data,
        )
    end

    return amplitude
end
