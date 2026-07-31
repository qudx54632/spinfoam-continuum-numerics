"""Single Euclidean EPRL 4-simplex vertex."""

const eprl_vertex_amplitude_cache = Dict()

empty!(eprl_vertex_amplitude_cache)

"""
Single Euclidean EPRL vertex in the same global tetrahedron basis as the BF
vertex.

The input `simplex` is an oriented 4-simplex.  The dictionaries `spin_data`
and `global_intertwiner_data` contain the diagonal SU(2) boundary labels.
"""
function eprl_vertex_amplitude(
    simplex,
    spin_data::AbstractDict,
    global_intertwiner_data::AbstractDict,
    gamma,
)
    0 <= gamma < 1 ||
        throw(ArgumentError("this first EPRL vertex implementation assumes 0 <= gamma < 1"))

    sigma = oriented_four_simplex(simplex)
    tetrahedra = oriented_simplex_tetrahedra(sigma)

    diagonal_spins = oriented_simplex_spins(sigma, spin_data)
    diagonal_intertwiners = Tuple(global_intertwiner_data[tau] for tau in tetrahedra)

    cache_key = (sigma, diagonal_spins, diagonal_intertwiners, gamma)
    if haskey(eprl_vertex_amplitude_cache, cache_key)
        return eprl_vertex_amplitude_cache[cache_key]
    end

    if !all(j -> eprl_spin_allowed(j, gamma), diagonal_spins)
        eprl_vertex_amplitude_cache[cache_key] = 0.0
        return 0.0
    end

    pair_ranges = [
        allowed_eprl_intertwiner_pairs(tau, spin_data, i, gamma)
        for (tau, i) in zip(tetrahedra, diagonal_intertwiners)
    ]

    if any(isempty, pair_ranges)
        eprl_vertex_amplitude_cache[cache_key] = 0.0
        return 0.0
    end

    left_spin_data = eprl_left_spin_data(spin_data, gamma)
    right_spin_data = eprl_right_spin_data(spin_data, gamma)

    total = 0.0

    for pair_tuple in Iterators.product(pair_ranges...)
        left_intertwiner_data = Dict()
        right_intertwiner_data = Dict()
        fusion_product = 1.0

        for a in 1:5
            tau = tetrahedra[a]
            i = diagonal_intertwiners[a]
            iL, iR = pair_tuple[a]

            left_intertwiner_data[tau] = iL
            right_intertwiner_data[tau] = iR

            spins = global_intertwiner_spins(tau, spin_data)
            fusion_product *= eprl_fusion_coefficient(spins, i, iL, iR, gamma)

            fusion_product == 0.0 && break
        end

        fusion_product == 0.0 && continue

        left_vertex = simplex_amplitude_recoupled_oriented_6j(
            sigma,
            left_spin_data,
            left_intertwiner_data,
        )

        left_vertex == 0.0 && continue

        right_vertex = simplex_amplitude_recoupled_oriented_6j(
            sigma,
            right_spin_data,
            right_intertwiner_data,
        )

        total += fusion_product * left_vertex * right_vertex
    end

    eprl_vertex_amplitude_cache[cache_key] = total

    return total
end
