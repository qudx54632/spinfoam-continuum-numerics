using WignerSymbols

const local_vertex_amplitude_cache = Dict()
const wigner3j_cache = Dict()
const magnetic_values_cache = Dict()

empty!(local_vertex_amplitude_cache)
empty!(wigner3j_cache)
empty!(magnetic_values_cache)

"""Allowed magnetic labels `m=-j,...,j`."""
function magnetic_values(j)
    if haskey(magnetic_values_cache, j)
        return magnetic_values_cache[j]
    end

    values = collect((-j):1:j)
    magnetic_values_cache[j] = values

    return values
end

"""Check whether `m` is an allowed magnetic label for spin `j`."""
allowed_magnetic_label(j, m) =
    -j <= m <= j && isinteger(j - m)

"""Cached Wigner 3j symbol, returning zero for inadmissible magnetic data."""
function cached_wigner3j(j1, j2, j3, m1, m2, m3)
    allowed_magnetic_label(j1, m1) || return 0.0
    allowed_magnetic_label(j2, m2) || return 0.0
    allowed_magnetic_label(j3, m3) || return 0.0
    iszero(m1 + m2 + m3) || return 0.0
    admissible_triangle((j1, j2, j3)) || return 0.0

    key = (j1, j2, j3, m1, m2, m3)

    if haskey(wigner3j_cache, key)
        return wigner3j_cache[key]
    end

    value = Float64(wigner3j(j1, j2, j3, m1, m2, m3))
    wigner3j_cache[key] = value

    return value
end

"""
Normalized four-valent intertwiner in the pairing

    (j1,j2) -> i -> (j3,j4).

The two Wigner 3j symbols impose

    n = -m1 - m2 = m3 + m4,

so there is no explicit n-sum.
"""
function four_valent_intertwiner(spins, i, magnetic_labels)
    spins = Tuple(spins)
    magnetic_labels = Tuple(magnetic_labels)
    j1, j2, j3, j4 = spins
    m1, m2, m3, m4 = magnetic_labels

    admissible_4valent(spins, i) || return 0.0

    n = -m1 - m2
    allowed_magnetic_label(i, n) || return 0.0
    iszero(-n + m3 + m4) || return 0.0

    first_symbol = cached_wigner3j(j1, j2, i, m1, m2, n)
    first_symbol == 0.0 && return 0.0

    second_symbol = cached_wigner3j(i, j3, j4, -n, m3, m4)
    second_symbol == 0.0 && return 0.0

    return sign_from_integer_exponent(i - n) *
           sqrt(Float64(dimension(i))) *
           first_symbol *
           second_symbol
end

"""
Local BF 4-simplex vertex as a normalized magnetic-index contraction.

The local spin order is

    J = (j12,j13,j14,j15,j23,j24,j25,j34,j35,j45),

and the local intertwiner order is

    k_tau = (i1,i2,i3,i4,i5).

This is the correctness-first definition of the local 15j amplitude:
five normalized four-valent intertwiners are glued with the SU(2) invariant
metric `(-1)^(j-m)`.
"""
function local_vertex_amplitude(J, k_tau)
    length(J) == 10 || throw(ArgumentError("J must have length 10"))
    length(k_tau) == 5 || throw(ArgumentError("k_tau must have length 5"))

    cache_key = (J, k_tau)
    if haskey(local_vertex_amplitude_cache, cache_key)
        return local_vertex_amplitude_cache[cache_key]
    end

    if !admissible_4simplex_intertwiners(J, k_tau)
        local_vertex_amplitude_cache[cache_key] = 0
        return 0
    end

    j12, j13, j14, j15, j23, j24, j25, j34, j35, j45 = J
    js = (j12, j13, j14, j15, j23, j24, j25, j34, j35, j45)
    node_spins = intertwiner_spin_lists(J)

    total = 0.0

    for m12 in magnetic_values(j12),
        m13 in magnetic_values(j13),
        m14 in magnetic_values(j14),
        m23 in magnetic_values(j23),
        m24 in magnetic_values(j24),
        m34 in magnetic_values(j34)

        m15 = -m12 - m13 - m14
        m25 = -m23 - m24 + m12
        m35 =  m13 + m23 - m34
        m45 =  m14 + m24 + m34

        allowed_magnetic_label(j15, m15) || continue
        allowed_magnetic_label(j25, m25) || continue
        allowed_magnetic_label(j35, m35) || continue
        allowed_magnetic_label(j45, m45) || continue

        all_magnetic_labels =
            (m12, m13, m14, m15, m23, m24, m25, m34, m35, m45)

        node_magnetic_labels = (
            (-m12, -m13, -m14, -m15),
            (-m23, -m24,  m12, -m25),
            ( m13,  m23, -m34, -m35),
            ( m14, -m45,  m24,  m34),
            ( m15,  m25,  m35,  m45),
        )

        metric_phase = prod(
            sign_from_integer_exponent(j - m)
            for (j, m) in zip(js, all_magnetic_labels);
            init = 1,
        )

        node_product = 1.0

        for a in 1:5
            node_product *= four_valent_intertwiner(
                node_spins[a],
                k_tau[a],
                node_magnetic_labels[a],
            )

            node_product == 0.0 && break
        end

        total += metric_phase * node_product
    end

    local_vertex_amplitude_cache[cache_key] = total

    return total
end
