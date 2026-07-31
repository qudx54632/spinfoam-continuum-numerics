const simplex_amplitude_recoupled_oriented_cache = Dict()

empty!(simplex_amplitude_recoupled_oriented_cache)

"""
Orientation-aware BF 4-simplex vertex in the fixed global tetrahedron basis.

This follows

    A^glob_sigma(J; i_tau)
      = sum_{k_tau} prod_tau R^sigma_tau(i_tau, k_tau)
        A^loc_sigma(J; k_tau).
"""
function simplex_amplitude_recoupled_oriented(
    simplex,
    spin_data::AbstractDict,
    global_intertwiner_data::AbstractDict,
)
    sigma = oriented_four_simplex(simplex)
    tetrahedra = oriented_simplex_tetrahedra(sigma)

    J = oriented_simplex_spins(sigma, spin_data)
    i_tau = Tuple(global_intertwiner_data[tau] for tau in tetrahedra)

    cache_key = (sigma, J, i_tau)
    if haskey(simplex_amplitude_recoupled_oriented_cache, cache_key)
        return simplex_amplitude_recoupled_oriented_cache[cache_key]
    end

    k_ranges = [
        allowed_oriented_local_intertwiners(sigma, tau, spin_data)
        for tau in tetrahedra
    ]

    if any(isempty, k_ranges)
        simplex_amplitude_recoupled_oriented_cache[cache_key] = 0
        return 0
    end

    amplitude = 0

    for k_tau in Iterators.product(k_ranges...)
        R = recoupling_product(sigma, tetrahedra, spin_data, i_tau, k_tau)

        if R != 0
            amplitude += R * local_vertex_amplitude(J, k_tau)
        end
    end

    simplex_amplitude_recoupled_oriented_cache[cache_key] = amplitude

    return amplitude
end

"""Fixed-label BF state-sum summand using recoupled oriented vertices."""
function state_sum_amplitude_recoupled_oriented(
    simplices,
    spin_data::AbstractDict,
    global_intertwiner_data::AbstractDict,
    internal_triangles,
)
    amplitude = triangle_amplitude_product(internal_triangles, spin_data)

    for simplex in simplices
        amplitude *= simplex_amplitude_recoupled_oriented(
            simplex,
            spin_data,
            global_intertwiner_data,
        )
    end

    return amplitude
end
