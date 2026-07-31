include(joinpath(@__DIR__, "..", "load_eprl.jl"))

"""
Semiclassical fusion-map check from arXiv:0809.3718.

For homogeneous spins `(j0,j0,j0,j0)`, compute

    g(iL,iR) = sum_i f^i_{iL iR}(j0) psi(i,j0),

where

    psi(i,j0) = exp[-sqrt(3)/2 (i-i0)^2/i0 + i*pi/2 (i-i0)]
    i0 = 2 j0 / sqrt(3).

The expected peak is near

    iL0 = (1-gamma) i0 / 2,
    iR0 = (1+gamma) i0 / 2,

for 0 <= gamma < 1.
"""

regular_tetrahedron_i0(j0) = 2 * Float64(j0) / sqrt(3)

function semiclassical_intertwiner_wavepacket(i, j0; i0 = regular_tetrahedron_i0(j0))
    delta_i = Float64(i) - Float64(i0)

    return exp(
        -sqrt(3) * delta_i^2 / (2 * Float64(i0)) +
        im * pi * delta_i / 2
    )
end

function homogeneous_allowed_intertwiners(j)
    spins = (j, j, j, j)

    return [
        i
        for i in half_integer_spins(2 * j)
        if admissible_4valent(spins, i)
    ]
end

function semiclassical_fusion_map_g(iL, iR, j0, gamma)
    spins = (j0, j0, j0, j0)
    total = 0.0 + 0.0im

    for i in homogeneous_allowed_intertwiners(j0)
        total += eprl_fusion_coefficient(spins, i, iL, iR, gamma) *
                 semiclassical_intertwiner_wavepacket(i, j0)
    end

    return total
end

function semiclassical_fusion_map_grid(j0, gamma)
    0 <= gamma < 1 ||
        throw(ArgumentError("this check assumes 0 <= gamma < 1"))
    eprl_spin_allowed(j0, gamma) ||
        throw(ArgumentError("j0 does not map to half-integer left/right spins"))

    jL = eprl_left_spin(j0, gamma)
    jR = eprl_right_spin(j0, gamma)
    iL_values = homogeneous_allowed_intertwiners(jL)
    iR_values = homogeneous_allowed_intertwiners(jR)

    return [
        (
            iL = iL,
            iR = iR,
            g = semiclassical_fusion_map_g(iL, iR, j0, gamma),
            abs_g = abs(semiclassical_fusion_map_g(iL, iR, j0, gamma)),
        )
        for iL in iL_values
        for iR in iR_values
    ]
end

function nearest_value(values, x)
    return values[argmin(abs.(Float64.(values) .- Float64(x)))]
end

function print_semiclassical_fusion_map_check(j0, gamma; top_count = 8)
    i0 = regular_tetrahedron_i0(j0)
    iL0 = (1 - Float64(gamma)) * i0 / 2
    iR0 = (1 + Float64(gamma)) * i0 / 2
    jL = eprl_left_spin(j0, gamma)
    jR = eprl_right_spin(j0, gamma)
    iL_values = homogeneous_allowed_intertwiners(jL)
    iR_values = homogeneous_allowed_intertwiners(jR)

    grid = semiclassical_fusion_map_grid(j0, gamma)
    sorted_grid = sort(grid, by = row -> row.abs_g, rev = true)
    peak = first(sorted_grid)

    println("Semiclassical EPRL fusion-map check")
    println("gamma                  = ", gamma)
    println("j0                     = ", j0)
    println("jL, jR                 = ", (jL, jR))
    println("i0 = 2 j0 / sqrt(3)    = ", i0)
    println("predicted iL0, iR0     = ", (iL0, iR0))
    println("nearest allowed peak   = ", (nearest_value(iL_values, iL0), nearest_value(iR_values, iR0)))
    println("exact |g| peak         = ", (peak.iL, peak.iR))
    println("peak |g|               = ", peak.abs_g)
    println("peak phase             = ", angle(peak.g))
    println()
    println("largest |g(iL,iR)| values:")

    for row in Iterators.take(sorted_grid, top_count)
        println(
            "  iL=", row.iL,
            "  iR=", row.iR,
            "  |g|=", row.abs_g,
            "  phase=", angle(row.g),
        )
    end
end

# Small default example.
#
# Increase `j0` to move closer to the large-spin regime.  The exact check is
# much cheaper than evaluating a full EPRL 4-simplex vertex.
gamma = 1 // 3
j0 = 6 // 1

print_semiclassical_fusion_map_check(j0, gamma)
