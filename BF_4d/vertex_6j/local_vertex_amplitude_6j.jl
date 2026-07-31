using WignerSymbols

"""Fast candidate for the local SU(2) BF 4-simplex vertex.

This file does not replace the magnetic-index definition in
`local_vertex_amplitude.jl`.  It gives a faster 6j expression for the same
local convention, and should be validated against the 3j contraction before
being used elsewhere.
"""

const local_vertex_amplitude_6j_cache = Dict()
const wigner6j_fast_cache = Dict()

empty!(local_vertex_amplitude_6j_cache)
empty!(wigner6j_fast_cache)

"""Cached Wigner 6j symbol, returning zero when the four triangles fail."""
function safe_wigner6j_fast(a, b, c, d, e, f)
    admissible_triangle((a, b, c)) || return 0.0
    admissible_triangle((a, e, f)) || return 0.0
    admissible_triangle((d, b, f)) || return 0.0
    admissible_triangle((d, e, c)) || return 0.0

    key = (a, b, c, d, e, f)
    if haskey(wigner6j_fast_cache, key)
        return wigner6j_fast_cache[key]
    end

    value = Float64(wigner6j(a, b, c, d, e, f))
    wigner6j_fast_cache[key] = value

    return value
end

"""
Convention-matching sign between the raw first-kind 5x6j expression and the
local magnetic-index convention used in `local_vertex_amplitude`.

The labels are ordered as

    J     = (j12,j13,j14,j15,j23,j24,j25,j34,j35,j45),
    k_tau = (i1,i2,i3,i4,i5).

This phase is intentionally kept in one visible place, because it is the
only convention-sensitive part of the fast local expression.
"""
function local_15j_6j_phase(J, k_tau)
    j12, j13, j14, j15, j23, j24, j25, j34, j35, j45 = J
    i1, i2, i3, i4, i5 = k_tau

    exponent =
        2j12 + j13 + 2j14 + j15 + j23 +
        3j24 + 3j34 + j35 + i1 + i3 + i5

    return sign_from_integer_exponent(exponent)
end

"""Raw first-kind 15j expression as a sum of five Wigner 6j symbols."""
function local_15j_6j_raw(J, k_tau)
    j12, j13, j14, j15, j23, j24, j25, j34, j35, j45 = J
    i1, i2, i3, i4, i5 = k_tau

    x_lower = maximum((
        abs(i1 - j25), abs(i5 - j14),
        abs(j14 - i5), abs(j35 - i4),
        abs(i4 - j35), abs(i3 - j24),
        abs(j24 - i3), abs(j13 - i2),
        abs(i2 - j13), abs(i1 - j25),
    ))

    x_upper = minimum((
        i1 + j25, i5 + j14,
        j14 + i5, j35 + i4,
        i4 + j35, i3 + j24,
        j24 + i3, j13 + i2,
        i2 + j13, i1 + j25,
    ))

    x_min_twice = ceil(Int, 2 * x_lower)
    x_max_twice = floor(Int, 2 * x_upper)
    x_min_twice > x_max_twice && return 0.0

    total = 0.0

    for x_twice in x_min_twice:x_max_twice
        x = x_twice // 2

        total += Float64(dimension(x)) *
            safe_wigner6j_fast(i1,  j25, x, i5,  j14, j15) *
            safe_wigner6j_fast(j14, i5,  x, j35, i4,  j45) *
            safe_wigner6j_fast(i4,  j35, x, i3,  j24, j34) *
            safe_wigner6j_fast(j24, i3,  x, j13, i2,  j23) *
            safe_wigner6j_fast(i2,  j13, x, i1,  j25, j12)
    end

    normalization = sqrt(Float64(prod(dimension(i) for i in k_tau; init = 1)))

    return normalization * total
end

"""Fast 6j candidate for the same local vertex as `local_vertex_amplitude`."""
function local_vertex_amplitude_6j(J, k_tau)
    length(J) == 10 || throw(ArgumentError("J must have length 10"))
    length(k_tau) == 5 || throw(ArgumentError("k_tau must have length 5"))

    J = Tuple(J)
    k_tau = Tuple(k_tau)
    cache_key = (J, k_tau)

    if haskey(local_vertex_amplitude_6j_cache, cache_key)
        return local_vertex_amplitude_6j_cache[cache_key]
    end

    if !admissible_4simplex_intertwiners(J, k_tau)
        local_vertex_amplitude_6j_cache[cache_key] = 0.0
        return 0.0
    end

    value = local_15j_6j_phase(J, k_tau) * local_15j_6j_raw(J, k_tau)
    local_vertex_amplitude_6j_cache[cache_key] = value

    return value
end
