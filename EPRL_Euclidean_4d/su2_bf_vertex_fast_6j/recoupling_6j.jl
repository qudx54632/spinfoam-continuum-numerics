"""Fast 6j expression for four-valent recoupling coefficients.

This file gives the same coefficient as `recoupling_R` in
`su2_bf_vertex/recoupling.jl`, but avoids the magnetic-label overlap sum.
The old 3j-overlap definition remains the reference convention.
"""

const recoupling_R_6j_cache = Dict()

empty!(recoupling_R_6j_cache)

same_unordered_pair(a, b, c, d) =
    (a == c && b == d) || (a == d && b == c)

function pairing_kind_and_canonical_order(local_to_global_position)
    p = local_to_global_position

    if same_unordered_pair(p[1], p[2], 1, 2) &&
       same_unordered_pair(p[3], p[4], 3, 4)
        return 1, (1, 2, 3, 4)
    elseif same_unordered_pair(p[1], p[2], 3, 4) &&
           same_unordered_pair(p[3], p[4], 1, 2)
        return 1, (1, 2, 3, 4)
    elseif same_unordered_pair(p[1], p[2], 1, 3) &&
           same_unordered_pair(p[3], p[4], 2, 4)
        return 2, (1, 3, 2, 4)
    elseif same_unordered_pair(p[1], p[2], 2, 4) &&
           same_unordered_pair(p[3], p[4], 1, 3)
        return 2, (1, 3, 2, 4)
    elseif same_unordered_pair(p[1], p[2], 1, 4) &&
           same_unordered_pair(p[3], p[4], 2, 3)
        return 3, (1, 4, 2, 3)
    elseif same_unordered_pair(p[1], p[2], 2, 3) &&
           same_unordered_pair(p[3], p[4], 1, 4)
        return 3, (1, 4, 2, 3)
    else
        throw(ArgumentError("invalid four-leg permutation"))
    end
end

"""
Phase relating an arbitrary ordered local basis to its canonical pairing.

The elementary signs are exactly the Condon--Shortley signs of the normalized
four-valent intertwiner:

    swap inside a pair:  (-1)^(j_a+j_b+k),
    swap the two pairs: (-1)^(2k).
"""
function local_order_to_canonical_pairing_phase(global_spins, local_to_global_position, k)
    pairing_kind, canonical_order =
        pairing_kind_and_canonical_order(local_to_global_position)

    p = local_to_global_position
    current = canonical_order
    phase = 1.0

    if same_unordered_pair(p[1], p[2], canonical_order[3], canonical_order[4])
        phase *= sign_from_integer_exponent(2k)
        current = (
            canonical_order[3],
            canonical_order[4],
            canonical_order[1],
            canonical_order[2],
        )
    end

    if current[1] != p[1]
        phase *= sign_from_integer_exponent(
            global_spins[current[1]] + global_spins[current[2]] + k,
        )
        current = (current[2], current[1], current[3], current[4])
    end

    if current[3] != p[3]
        phase *= sign_from_integer_exponent(
            global_spins[current[3]] + global_spins[current[4]] + k,
        )
        current = (current[1], current[2], current[4], current[3])
    end

    current == p ||
        throw(ArgumentError("could not reduce local order to canonical pairing"))

    return pairing_kind, phase
end

"""Racah move from `(12)(34)` to `(13)(24)`."""
function recoupling_13_6j(global_spins, i, k)
    j1, j2, j3, j4 = global_spins

    return sign_from_integer_exponent(j2 + j3 + i + k) *
           sqrt(Float64(dimension(i) * dimension(k))) *
           safe_wigner6j_fast(j1, j2, i, j4, j3, k)
end

"""Racah move from `(12)(34)` to `(14)(23)`."""
function recoupling_14_6j(global_spins, i, k)
    j1, j2, j3, j4 = global_spins

    return sign_from_integer_exponent(i + j3 + j4) *
           recoupling_13_6j((j1, j2, j4, j3), i, k)
end

function canonical_recoupling_6j(global_spins, i, k, pairing_kind)
    if pairing_kind == 1
        return i == k ? 1.0 : 0.0
    elseif pairing_kind == 2
        return recoupling_13_6j(global_spins, i, k)
    elseif pairing_kind == 3
        return recoupling_14_6j(global_spins, i, k)
    else
        throw(ArgumentError("unknown pairing kind"))
    end
end

"""
Fast recoupling coefficient

    R^sigma_tau(i_tau,k_tau)
      = <i_tau^global | k_tau^local,sigma>.
"""
function recoupling_R_6j(simplex, tetrahedron, spin_data::AbstractDict, i_tau, k_tau)
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
    if haskey(recoupling_R_6j_cache, cache_key)
        return recoupling_R_6j_cache[cache_key]
    end

    local_to_global_position =
        Tuple(findfirst(==(key), global_keys) for key in local_keys)

    pairing_kind, order_phase = local_order_to_canonical_pairing_phase(
        global_spins,
        local_to_global_position,
        k_tau,
    )

    value =
        order_phase *
        canonical_recoupling_6j(global_spins, i_tau, k_tau, pairing_kind)

    recoupling_R_6j_cache[cache_key] = value

    return value
end

"""Product of the five fast recoupling coefficients in one 4-simplex."""
function recoupling_product_6j(simplex, tetrahedra, spin_data::AbstractDict, global_i, local_k)
    R = 1.0

    for a in 1:5
        R *= recoupling_R_6j(
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
