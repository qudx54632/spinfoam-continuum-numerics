const recoupling_R_cache = Dict()

empty!(recoupling_R_cache)

"""Return the local spin order attached to one tetrahedron of a simplex."""
function oriented_local_intertwiner_spin_keys(simplex, tetrahedron)
    tetrahedra = oriented_simplex_intertwiner_keys(simplex)
    position = findfirst(==(tetrahedron_key(tetrahedron...)), tetrahedra)

    position === nothing &&
        throw(ArgumentError("the tetrahedron is not in this oriented 4-simplex"))

    local_spin_keys = oriented_simplex_spin_keys(simplex)

    return intertwiner_spin_lists(local_spin_keys)[position]
end

"""Allowed local channel labels for one tetrahedron of a simplex."""
function allowed_oriented_local_intertwiners(simplex, tetrahedron, spin_data::AbstractDict)
    spin_keys = oriented_local_intertwiner_spin_keys(simplex, tetrahedron)
    spins = Tuple(spin_data[key] for key in spin_keys)
    j1, j2, j3, j4 = spins
    i_max = min(j1 + j2, j3 + j4)

    return [
        k
        for k in half_integer_spins(i_max)
        if admissible_4valent(spins, k)
    ]
end

"""
Recoupling coefficient

    R^sigma_tau(i_tau, k_tau)
      = <i_tau^global | k_tau^local,sigma>.

The ordered global triangle list is

    global_keys = (g1,g2,g3,g4).

The ordered local triangle list is the same four triangles, possibly in a
different order,

    local_keys = (g_{p1},g_{p2},g_{p3},g_{p4}).

The code below is exactly the overlap formula in the note: sum over magnetic
labels in the global order, then reorder them by `(p1,p2,p3,p4)` before
evaluating the local intertwiner.
"""
function recoupling_R(simplex, tetrahedron, spin_data::AbstractDict, i_tau, k_tau)
    global_keys = global_intertwiner_spin_keys(tetrahedron)
    local_keys = oriented_local_intertwiner_spin_keys(simplex, tetrahedron)

    Set(global_keys) == Set(local_keys) ||
        throw(ArgumentError("global and local bases must use the same four triangles"))

    global_spins = Tuple(spin_data[key] for key in global_keys)
    local_spins = Tuple(spin_data[key] for key in local_keys)

    if !admissible_4valent(global_spins, i_tau) ||
       !admissible_4valent(local_spins, k_tau)
        return 0.0
    end

    cache_key = (global_keys, local_keys, global_spins, i_tau, k_tau)
    if haskey(recoupling_R_cache, cache_key)
        return recoupling_R_cache[cache_key]
    end

    if global_keys == local_keys
        value = i_tau == k_tau ? 1.0 : 0.0
        recoupling_R_cache[cache_key] = value

        return value
    end

    local_to_global_position = Tuple(findfirst(==(key), global_keys) for key in local_keys)
    j1, j2, j3, j4 = global_spins
    total = 0.0

    for m1 in magnetic_values(j1),
        m2 in magnetic_values(j2),
        m3 in magnetic_values(j3)

        m4 = -m1 - m2 - m3
        allowed_magnetic_label(j4, m4) || continue

        global_magnetic_labels = (m1, m2, m3, m4)
        local_magnetic_labels =
            Tuple(global_magnetic_labels[position] for position in local_to_global_position)

        global_intertwiner = four_valent_intertwiner(
            global_spins,
            i_tau,
            global_magnetic_labels,
        )

        local_intertwiner = four_valent_intertwiner(
            local_spins,
            k_tau,
            local_magnetic_labels,
        )

        total += conj(global_intertwiner) * local_intertwiner
    end

    recoupling_R_cache[cache_key] = total

    return total
end

"""Product of the five recoupling coefficients in one 4-simplex."""
function recoupling_product(simplex, tetrahedra, spin_data::AbstractDict, global_i, local_k)
    R = 1.0

    for a in 1:5
        R *= recoupling_R(
            simplex,
            tetrahedra[a],
            spin_data,
            global_i[a],
            local_k[a],
        )

        R == 0.0 && return 0.0
    end

    return R
end
