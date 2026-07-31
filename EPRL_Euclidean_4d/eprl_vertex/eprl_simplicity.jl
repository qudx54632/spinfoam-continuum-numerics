"""Simplicity map for the Euclidean EPRL vertex."""

"""Left SU(2) spin `jL = |1 - gamma| j / 2`."""
eprl_left_spin(j, gamma) = abs(1 - gamma) * j / 2

"""Right SU(2) spin `jR = (1 + gamma) j / 2`."""
eprl_right_spin(j, gamma) = (1 + gamma) * j / 2

"""Check whether `j` is a non-negative half-integer."""
is_half_integer_spin(j) = j >= 0 && isinteger(2 * j)

"""Check whether one diagonal spin has admissible left and right images."""
function eprl_spin_allowed(j, gamma)
    jL = eprl_left_spin(j, gamma)
    jR = eprl_right_spin(j, gamma)

    return is_half_integer_spin(j) &&
           is_half_integer_spin(jL) &&
           is_half_integer_spin(jR)
end

"""Map every triangle spin in `spin_data` to the left SU(2) sector."""
function eprl_left_spin_data(spin_data::AbstractDict, gamma)
    return Dict(key => eprl_left_spin(j, gamma) for (key, j) in spin_data)
end

"""Map every triangle spin in `spin_data` to the right SU(2) sector."""
function eprl_right_spin_data(spin_data::AbstractDict, gamma)
    return Dict(key => eprl_right_spin(j, gamma) for (key, j) in spin_data)
end

