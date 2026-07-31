include(joinpath(@__DIR__, "semiclassical", "semiclassical_vertex_tools.jl"))

using Plots

# Single-vertex coherent-boundary Regge check.
#
# Change the boundary data directly in this file.

simplex = (1, 2, 3, 4, 5)
gamma = 1 // 3
base_spin = 3 // 2

# This scans lambda = 1, ..., lambda_stop.  The coherent vertex sums all five
# boundary intertwiners and gets expensive quickly, so reduce this range while
# debugging.
lambda_stop = parse(Int, get(ENV, "EPRL_LAMBDA_STOP", "3"))
lambda_values = 1:lambda_stop

# Default: regular 4-simplex boundary data with all ten spins equal to `j`.
#
# To test another boundary state, replace this function by one returning
#
#     spin_data, normal_data
#
# where `spin_data` is keyed by triangle and `normal_data` is keyed by
# `(tetrahedron, triangle)`.
boundary_data(simplex, j) = regular_4simplex_boundary_data(simplex, j)

function save_semiclassical_scan_csv(filename, scan)
    open(filename, "w") do io
        println(
            io,
            "lambda,j,abs_W,arg_W,real_W,imag_W,real_lambda12_W,imag_lambda12_W,S_int,S_ext,closure_error",
        )

        for row in scan
            println(
                io,
                row.lambda, ",",
                Float64(row.j), ",",
                abs(row.W), ",",
                angle(row.W), ",",
                real(row.W), ",",
                imag(row.W), ",",
                real(row.scaled_W), ",",
                imag(row.scaled_W), ",",
                row.S_int, ",",
                row.S_ext, ",",
                row.closure_error,
            )
        end
    end
end

scan = []

println("Euclidean EPRL coherent single-vertex semiclassical check")
println("simplex                = ", simplex)
println("gamma                  = ", gamma)
println("base spin              = ", base_spin)
println("lambda values          = ", collect(lambda_values))
println()
println("lambda  j       |W|                 arg(W)       Re(W)                 Im(W)")

for lambda in lambda_values
    j = lambda * base_spin
    spin_data, normal_data = boundary_data(simplex, j)

    closure_error = max_closure_norm(simplex, spin_data, normal_data)
    W = coherent_single_vertex_amplitude(simplex, spin_data, normal_data, gamma)
    scaled_W = lambda^12 * W

    S_int = regular_4simplex_regge_action(j, gamma; dihedral = :interior)
    S_ext = regular_4simplex_regge_action(j, gamma; dihedral = :exterior)

    push!(
        scan,
        (
            lambda = lambda,
            j = j,
            W = W,
            scaled_W = scaled_W,
            S_int = S_int,
            S_ext = S_ext,
            closure_error = closure_error,
        ),
    )

    println(
        lambda, "       ",
        Float64(j), "     ",
        abs(W), "    ",
        angle(W), "    ",
        real(W), "    ",
        imag(W),
    )
end

output_dir = joinpath(@__DIR__, "outputs")
mkpath(output_dir)

csv_file = joinpath(output_dir, "semiclassical_single_vertex_scan.csv")
plot_file = joinpath(output_dir, "semiclassical_single_vertex_W.png")

save_semiclassical_scan_csv(csv_file, scan)

x = [row.lambda for row in scan]
real_W = [real(row.W) for row in scan]
imag_W = [imag(row.W) for row in scan]

p = plot(
    x,
    real_W;
    marker = :circle,
    label = "Re(W)",
    xlabel = "lambda",
    ylabel = "W",
)

plot!(
    p,
    x,
    imag_W;
    marker = :diamond,
    label = "Im(W)",
)

savefig(p, plot_file)

println()
println("Expected large-spin structure:")
println("  W(lambda) ~ lambda^(-12) [N_+ exp(i S_Regge) + N_- exp(-i S_Regge)]")
println("Therefore lambda^12 W(lambda) should show Regge-frequency oscillations.")
println()
println("saved data             = ", csv_file)
println("saved plot             = ", plot_file)
