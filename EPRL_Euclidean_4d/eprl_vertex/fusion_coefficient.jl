"""Fusion coefficients for the Euclidean EPRL vertex.

This file implements the `0 <= gamma < 1` case, where

    j = jL + jR.

The 9j-symbol entering the fusion coefficient is evaluated using the
special formula quoted in Appendix A of arXiv:0809.3718.
"""

const eprl_q_cache = Dict()
const eprl_fusion_cache = Dict()

empty!(eprl_q_cache)
empty!(eprl_fusion_cache)

function nonnegative_integer_value(x)
    isinteger(x) || return nothing
    n = Int(x)
    n >= 0 || return nothing

    return n
end

function log_factorial_integer(n::Integer)
    n >= 0 || throw(ArgumentError("factorial argument must be non-negative"))

    total = 0.0
    for k in 2:n
        total += log(k)
    end

    return total
end

function log_factorial_spin_argument(x)
    n = nonnegative_integer_value(x)
    n === nothing && return nothing

    return log_factorial_integer(n)
end

function log_factorial_product(args)
    total = 0.0

    for arg in args
        value = log_factorial_spin_argument(arg)
        value === nothing && return nothing
        total += value
    end

    return total
end

"""
Special 9j-symbol

    { a  f  c }
    { b  g  d }
    { a+b h c+d }

with the convention of Appendix A, Eq. (37), of arXiv:0809.3718.
"""
function special_wigner9j_sum_columns(a, f, c, b, g, d, h)
    admissible_triangle((a, c, f)) || return 0.0
    admissible_triangle((b, d, g)) || return 0.0
    admissible_triangle((a + b, c + d, h)) || return 0.0

    m_f = a - c
    m_g = b - d
    m_h = -(a + b - (c + d))

    three_j = cached_wigner3j(f, g, h, m_f, m_g, m_h)
    three_j == 0.0 && return 0.0

    numerator_args = (
        2a,
        2b,
        2c,
        2d,
        a + b + c + d - h,
        a + b + c + d + h + 1,
    )

    denominator_args = (
        2a + 2b + 1,
        2c + 2d + 1,
        a + c - f,
        a + c + f + 1,
        b + d - g,
        b + d + g + 1,
    )

    log_num = log_factorial_product(numerator_args)
    log_den = log_factorial_product(denominator_args)

    (log_num === nothing || log_den === nothing) && return 0.0

    phase = sign_from_integer_exponent(f - g + a + b - (c + d))
    prefactor = exp(0.5 * (log_num - log_den))

    return phase * three_j * prefactor
end

"""
The q-coefficient

    q^i_{iL iR}(ja,jb; gamma)
      =
        { jaL  iL  jbL }
        { jaR  iR  jbR }
        { ja   i   jb  }.

Only `0 <= gamma < 1` is implemented in this first version.
"""
function eprl_q_coefficient(ja, jb, i, iL, iR, gamma)
    0 <= gamma < 1 ||
        throw(ArgumentError("this q-coefficient currently assumes 0 <= gamma < 1"))

    key = (ja, jb, i, iL, iR, gamma)
    if haskey(eprl_q_cache, key)
        return eprl_q_cache[key]
    end

    jaL = eprl_left_spin(ja, gamma)
    jaR = eprl_right_spin(ja, gamma)
    jbL = eprl_left_spin(jb, gamma)
    jbR = eprl_right_spin(jb, gamma)

    if !all(is_half_integer_spin, (jaL, jaR, jbL, jbR))
        eprl_q_cache[key] = 0.0
        return 0.0
    end

    value = special_wigner9j_sum_columns(jaL, iL, jbL, jaR, iR, jbR, i)
    eprl_q_cache[key] = value

    return value
end

"""
Fusion coefficient in the fixed global tetrahedron basis:

    f^i_{iL iR}(j1,j2,j3,j4; gamma).
"""
function eprl_fusion_coefficient(spins, i, iL, iR, gamma)
    spins = Tuple(spins)
    length(spins) == 4 || throw(ArgumentError("spins must have length 4"))
    j1, j2, j3, j4 = spins

    key = (spins, i, iL, iR, gamma)
    if haskey(eprl_fusion_cache, key)
        return eprl_fusion_cache[key]
    end

    if !admissible_4valent(spins, i) ||
       !admissible_triangle((iL, iR, i))
        eprl_fusion_cache[key] = 0.0
        return 0.0
    end

    left_spins = Tuple(eprl_left_spin(j, gamma) for j in spins)
    right_spins = Tuple(eprl_right_spin(j, gamma) for j in spins)

    if !all(is_half_integer_spin, left_spins) ||
       !all(is_half_integer_spin, right_spins) ||
       !admissible_4valent(left_spins, iL) ||
       !admissible_4valent(right_spins, iR)
        eprl_fusion_cache[key] = 0.0
        return 0.0
    end

    q12 = eprl_q_coefficient(j1, j2, i, iL, iR, gamma)
    q34 = eprl_q_coefficient(j3, j4, i, iL, iR, gamma)

    if q12 == 0.0 || q34 == 0.0
        eprl_fusion_cache[key] = 0.0
        return 0.0
    end

    phase = sign_from_integer_exponent(j1 - j2 + j3 - j4)
    normalization = sqrt(Float64(dimension(i) * dimension(iL) * dimension(iR) *
                                 prod(dimension(j) for j in spins)))

    value = phase * normalization * q12 * q34
    eprl_fusion_cache[key] = value

    return value
end

"""Allowed `(iL,iR)` pairs for one tetrahedron and one diagonal intertwiner."""
function allowed_eprl_intertwiner_pairs(tetrahedron, spin_data::AbstractDict, i, gamma)
    spins = global_intertwiner_spins(tetrahedron, spin_data)
    left_spins = Tuple(eprl_left_spin(j, gamma) for j in spins)
    right_spins = Tuple(eprl_right_spin(j, gamma) for j in spins)

    if !admissible_4valent(spins, i) ||
       !all(is_half_integer_spin, left_spins) ||
       !all(is_half_integer_spin, right_spins)
        return Tuple{Real,Real}[]
    end

    j1L, j2L, j3L, j4L = left_spins
    j1R, j2R, j3R, j4R = right_spins
    iL_max = min(j1L + j2L, j3L + j4L)
    iR_max = min(j1R + j2R, j3R + j4R)

    pairs = []

    for iL in half_integer_spins(iL_max), iR in half_integer_spins(iR_max)
        if admissible_4valent(left_spins, iL) &&
           admissible_4valent(right_spins, iR) &&
           admissible_triangle((iL, iR, i))
            push!(pairs, (iL, iR))
        end
    end

    return pairs
end

