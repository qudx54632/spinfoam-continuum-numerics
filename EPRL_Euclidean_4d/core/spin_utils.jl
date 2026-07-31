"""Basic SU(2) spin utilities used by the oriented 4D BF state sum."""

"""SU(2) representation dimension `d_j = 2j + 1`."""
dimension(j::Real) = 2 * j + 1

"""Return all half-integer spins from `0` through `j_max`."""
function half_integer_spins(j_max::Real)
    j_max >= 0 || throw(ArgumentError("j_max must be nonnegative"))
    isinteger(2 * j_max) || throw(ArgumentError("j_max must be a half-integer"))

    return [n // 2 for n in 0:Int(2 * j_max)]
end

"""Test the SU(2) admissibility conditions for three spins."""
function admissible_triangle(spins)
    j1, j2, j3 = spins

    return all(j -> j >= 0 && isinteger(2 * j), spins) &&
           j1 + j2 >= j3 &&
           j1 + j3 >= j2 &&
           j2 + j3 >= j1 &&
           isinteger(j1 + j2 + j3)
end

"""
Test a four-valent intertwiner channel.

The convention is

    (j1, j2) -> i -> (j3, j4).
"""
function admissible_4valent(spins, i)
    j1, j2, j3, j4 = spins

    return admissible_triangle((j1, j2, i)) &&
           admissible_triangle((i, j3, j4))
end

"""
Return the five four-valent spin lists in the local magnetic-index 15j order.

The input order is

    (j12, j13, j14, j15, j23, j24, j25, j34, j35, j45).
"""
function intertwiner_spin_lists(js)
    j12, j13, j14, j15, j23, j24, j25, j34, j35, j45 = js

    return (
        (j12, j13, j14, j15),
        (j23, j24, j12, j25),
        (j13, j23, j34, j35),
        (j14, j45, j24, j34),
        (j15, j25, j35, j45),
    )
end

"""Test whether the five intertwiners of a local 4-simplex are admissible."""
function admissible_4simplex_intertwiners(js, intertwiners)
    length(js) == 10 || throw(ArgumentError("js must have length 10"))
    length(intertwiners) == 5 || throw(ArgumentError("intertwiners must have length 5"))

    return all(
        pair -> admissible_4valent(pair[1], pair[2]),
        zip(intertwiner_spin_lists(js), intertwiners),
    )
end

"""Return `(-1)^n`, checking that `n` is an integer."""
function sign_from_integer_exponent(n)
    isinteger(n) || throw(ArgumentError("the phase exponent must be an integer"))

    return (-1)^Int(n)
end
