using Distributed

const eprl_dir = @__DIR__

@everywhere begin
    include(joinpath($eprl_dir, "semiclassical", "semiclassical_vertex_tools.jl"))
end

using Plots

const simplex = (1, 2, 3, 4, 5)
const base_spin = 1 // 2

# Equispin scaling:
#
#     j_ab(lambda) = lambda/2.
#
# This is the BF vertex amplitude contracted with Livine-Speziale coherent
# boundary intertwiners.  It has no internal-spin cutoff; the only scan
# parameter here is lambda.
lambda_start = parse(Int, get(ENV, "SU2_BF_LAMBDA_START", "1"))
lambda_stop = parse(Int, get(ENV, "SU2_BF_LAMBDA_STOP", "70"))
lambda_values = collect(lambda_start:lambda_stop)

@everywhere function su2_bf_coherent_lambda_row(lambda::Integer, base_spin, simplex)
    j = lambda * base_spin
    spin_data, normal_data = regular_4simplex_boundary_data(simplex, j)

    closure_error = max_closure_norm(simplex, spin_data, normal_data)

    elapsed = @elapsed begin
        W = coherent_su2_bf_single_vertex_amplitude_6j(
            simplex,
            spin_data,
            normal_data,
        )
    end

    scaled_W = lambda^6 * W
    S_int = regular_4simplex_su2_bf_regge_action(j; dihedral = :interior)
    S_ext = regular_4simplex_su2_bf_regge_action(j; dihedral = :exterior)

    return (
        lambda = lambda,
        j = j,
        W = W,
        scaled_W = scaled_W,
        S_int = S_int,
        S_ext = S_ext,
        closure_error = closure_error,
        elapsed = elapsed,
    )
end

function warmup_su2_bf_workers()
    if nprocs() == 1
        su2_bf_coherent_lambda_row(1, base_spin, simplex)
    else
        @sync for worker in workers()
            @async remotecall_fetch(
                su2_bf_coherent_lambda_row,
                worker,
                1,
                base_spin,
                simplex,
            )
        end
    end

    return nothing
end

function run_su2_bf_scan(lambda_values)
    if nprocs() == 1
        return [su2_bf_coherent_lambda_row(lambda, base_spin, simplex)
                for lambda in lambda_values]
    end

    return pmap(
        lambda -> su2_bf_coherent_lambda_row(lambda, base_spin, simplex),
        lambda_values,
    )
end

function save_su2_bf_scan_csv(filename, scan)
    open(filename, "w") do io
        println(
            io,
            "lambda,j,abs_W,arg_W,real_W,imag_W,real_lambda6_W,imag_lambda6_W,abs_lambda6_W,S_int,S_ext,closure_error,elapsed_seconds",
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
                abs(row.scaled_W), ",",
                row.S_int, ",",
                row.S_ext, ",",
                row.closure_error, ",",
                row.elapsed,
            )
        end
    end
end

function plot_su2_bf_scan(scan, plot_file)
    x = [row.lambda for row in scan]
    scaled_real = [real(row.scaled_W) for row in scan]
    scaled_imag = [imag(row.scaled_W) for row in scan]
    scaled_abs = [abs(row.scaled_W) for row in scan]
    cos_regge = [cos(row.S_ext) for row in scan]

    scale = maximum(abs.(vcat(scaled_real, scaled_imag)))
    scaled_cos_regge = scale == 0 ? cos_regge : scale .* cos_regge

    p1 = plot(
        x,
        scaled_real;
        marker = :circle,
        label = "Re(lambda^6 W)",
        ylabel = "lambda^6 W",
    )

    plot!(
        p1,
        x,
        scaled_imag;
        marker = :diamond,
        label = "Im(lambda^6 W)",
    )

    plot!(
        p1,
        x,
        scaled_cos_regge;
        linestyle = :dash,
        label = "rescaled cos(S_BF exterior)",
    )

    p2 = plot(
        x,
        scaled_abs;
        marker = :circle,
        label = "|lambda^6 W|",
        xlabel = "lambda",
        ylabel = "|lambda^6 W|",
    )

    plot(p1, p2; layout = (2, 1), size = (900, 700))
    savefig(plot_file)
end

println("SU(2) BF vertex amplitude with coherent boundary state")
println("vertex backend = fast 6j local vertex + fast 6j recoupling")
println("simplex        = ", simplex)
println("base spin      = ", base_spin)
println("lambda values  = ", lambda_values)
println("processes      = ", nprocs(), " total, ", nworkers(), " workers")
println("cutoff         = none; this is a fixed-boundary coherent vertex amplitude")
println()

println("Warming up workers...")
warmup_su2_bf_workers()

println("Running lambda scan...")
scan = sort(run_su2_bf_scan(lambda_values); by = row -> row.lambda)

println()
println("lambda  j       |W|                 arg(W)       Re(lambda^6 W)       Im(lambda^6 W)       time")
for row in scan
    println(
        row.lambda, "       ",
        Float64(row.j), "     ",
        abs(row.W), "    ",
        angle(row.W), "    ",
        real(row.scaled_W), "    ",
        imag(row.scaled_W), "    ",
        row.elapsed,
    )
end

output_dir = joinpath(@__DIR__, "outputs")
mkpath(output_dir)

csv_file = joinpath(output_dir, "su2_bf_coherent_vertex_lambda_scan.csv")
plot_file = joinpath(output_dir, "su2_bf_coherent_lambda6_scan.png")

save_su2_bf_scan_csv(csv_file, scan)
plot_su2_bf_scan(scan, plot_file)

println()
println("Expected SU(2) BF large-spin structure:")
println("  W(lambda) ~ lambda^(-6) cos(lambda S_BF + constant phase)")
println("To see the cosine sign oscillation, use Re(lambda^6 W) or Im(lambda^6 W).")
println("|lambda^6 W| is useful as a magnitude/envelope check, but it removes signs.")
println()
println("saved data     = ", csv_file)
println("saved plot     = ", plot_file)
