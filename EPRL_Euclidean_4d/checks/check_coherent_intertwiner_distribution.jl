include(joinpath(@__DIR__, "..", "load_eprl.jl"))

using LinearAlgebra
using Plots
using Printf

# Distribution of one Livine-Speziale coherent intertwiner.
#
# This is a diagnostic plot for comparing with Fig. 7 of arXiv:1903.12624.
# It is not used by the vertex-amplitude calculation.  The amplitude code
# keeps using the normalized-basis coefficient returned by
# coherent_intertwiner_coefficient.
#
# This reproduces the standard elementary check:
#
#   j_ab = 10,
#   normals = normals of an equilateral tetrahedron,
#   coefficient c_i as a function of i = j12.
#
# For four equal spins j, the admissible channel is i = 0, 1, ..., 2j.

j = parse(Rational{Int}, get(ENV, "COHERENT_INTERTWINER_J", "10//1"))
spins = (j, j, j, j)

# Outward unit normals of a regular tetrahedron in R^3.
r = sqrt(3.0)
normals = (
    ( 1.0 / r,  1.0 / r,  1.0 / r),
    ( 1.0 / r, -1.0 / r, -1.0 / r),
    (-1.0 / r,  1.0 / r, -1.0 / r),
    (-1.0 / r, -1.0 / r,  1.0 / r),
)

i_values = [
    i
    for i in half_integer_spins(2j)
    if admissible_4valent(spins, i)
]
normalized_basis_coefficients = [
    coherent_intertwiner_coefficient(spins, i, normals)
    for i in i_values
]

# The code internally uses normalized four-valent intertwiners.  The paper's
# Fig. 7 uses the coefficient c_i defined with the unnormalized four-legged
# Wigner symbol.  These two conventions differ by sqrt(d_i):
#
#     c_i^(normalized basis) = sqrt(2i+1) c_i^(paper).
paper_coefficients = [
    abs(c) / sqrt(Float64(dimension(i)))
    for (i, c) in zip(i_values, normalized_basis_coefficients)
]

normalized_basis_weights = abs2.(normalized_basis_coefficients)
normalized_basis_probability =
    normalized_basis_weights ./ sum(normalized_basis_weights)

# The classical value of the recoupled vector |J_1 + J_2|.
i_classical = Float64(j) *
    norm(Float64.(normals[1]) .+ Float64.(normals[2]))

peak_index = argmax(paper_coefficients)
peak_i = i_values[peak_index]
peak_coefficient = paper_coefficients[peak_index]

output_dir = joinpath(@__DIR__, "..", "outputs", "checks")
mkpath(output_dir)

plot_file = joinpath(output_dir, "coherent_intertwiner_j$(Int(j))_distribution.png")
csv_file = joinpath(output_dir, "coherent_intertwiner_j$(Int(j))_distribution.csv")

open(csv_file, "w") do io
    println(
        io,
        "i,real_c_normalized_basis,imag_c_normalized_basis,abs_c_normalized_basis,c_paper,abs2_c_normalized_basis,normalized_abs2_c_normalized_basis",
    )

    for (i, c, cp, w, nw) in zip(
        i_values,
        normalized_basis_coefficients,
        paper_coefficients,
        normalized_basis_weights,
        normalized_basis_probability,
    )
        println(
            io,
            Float64(i), ",",
            real(c), ",",
            imag(c), ",",
            abs(c), ",",
            cp, ",",
            w, ",",
            nw,
        )
    end
end

p = plot(
    Float64.(i_values),
    paper_coefficients;
    seriestype = :bar,
    label = "paper convention c_i",
    xlabel = "intertwiner label i = j12",
    ylabel = "c_i",
    title = "Equilateral coherent intertwiner, j_ab = $j",
    size = (850, 550),
    legend = :topright,
)

vline!(
    p,
    [i_classical];
    linestyle = :dash,
    linewidth = 2,
    label = "classical |J1 + J2|",
)

scatter!(
    p,
    [Float64(peak_i)],
    [peak_coefficient];
    marker = :circle,
    markersize = 6,
    label = "peak",
)

savefig(p, plot_file)

println("Coherent intertwiner distribution check")
println("spins       = ", spins)
println("i range     = ", first(i_values), ":", last(i_values))
@printf("classical i = %.12f\n", i_classical)
println("peak i      = ", peak_i)
@printf("peak c_i    = %.12g\n", peak_coefficient)
println("saved CSV   = ", csv_file)
println("saved plot  = ", plot_file)
