"""Boundary spin-network phase used when comparing before and after sides."""
function pachner_boundary_orientation_phase(boundary_spins)
    return sign_from_integer_exponent(2 * sum(boundary_spins))
end
