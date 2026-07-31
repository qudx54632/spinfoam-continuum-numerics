using Distributed

const eprl_dir = @__DIR__

function add_requested_su2_bf_workers()
    requested_workers = parse(Int, get(ENV, "SU2_BF_NWORKERS", "0"))
    requested_workers <= 0 && return nothing

    if nprocs() > 1
        @warn(
            "Workers are already running. SU2_BF_NWORKERS is ignored. " *
            "For large runs, prefer `SU2_BF_NWORKERS=N julia main_check_su2_bf_semiclassical.jl` " *
            "instead of `julia -p N ...`, so the script can use topology=:master_worker."
        )
        return nothing
    end

    active_project = Base.active_project()

    if active_project === nothing
        addprocs(requested_workers; topology = :master_worker)
    else
        addprocs(
            requested_workers;
            topology = :master_worker,
            exeflags = "--project=$(active_project)",
        )
    end

    return nothing
end

add_requested_su2_bf_workers()

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
chunks_per_worker = parse(Int, get(ENV, "SU2_BF_CHUNKS_PER_WORKER", "2"))

@everywhere function su2_bf_intertwiner_chunk_sum(
    simplex,
    spin_data::AbstractDict,
    tetrahedra,
    coefficient_table::AbstractDict,
    intertwiner_ranges,
    first_two_intertwiner_pairs,
)
    sigma = oriented_four_simplex(simplex)
    total = 0.0 + 0.0im

    for (i1, i2) in first_two_intertwiner_pairs,
        i3 in intertwiner_ranges[3],
        i4 in intertwiner_ranges[4],
        i5 in intertwiner_ranges[5]

        intertwiners = (i1, i2, i3, i4, i5)
        coefficient_product = coherent_intertwiner_coefficient_product(
            coefficient_table,
            tetrahedra,
            intertwiners,
        )

        abs(coefficient_product) == 0.0 && continue

        intertwiner_data = Dict(
            tetrahedron => i
            for (tetrahedron, i) in zip(tetrahedra, intertwiners)
        )

        total += simplex_amplitude_recoupled_oriented_6j(
            sigma,
            spin_data,
            intertwiner_data,
        ) * coefficient_product
    end

    return total
end

@everywhere function warmup_su2_bf_worker(base_spin, simplex)
    j = base_spin
    spin_data, normal_data = regular_4simplex_boundary_data(simplex, j)

    return coherent_su2_bf_single_vertex_amplitude_6j(
        simplex,
        spin_data,
        normal_data,
    )
end

function split_into_chunks(values, number_of_chunks)
    number_of_chunks = max(1, min(number_of_chunks, length(values)))
    chunk_size = cld(length(values), number_of_chunks)

    return [
        values[first:min(first + chunk_size - 1, length(values))]
        for first in 1:chunk_size:length(values)
    ]
end

function su2_bf_coherent_amplitude_parallel_over_intertwiners(
    simplex,
    spin_data::AbstractDict,
    normal_data::AbstractDict,
)
    sigma = oriented_four_simplex(simplex)
    tetrahedra = oriented_simplex_tetrahedra(sigma)

    intertwiner_ranges = [
        allowed_global_intertwiners(tetrahedron, spin_data)
        for tetrahedron in tetrahedra
    ]

    any(isempty, intertwiner_ranges) && return 0.0 + 0.0im

    coefficient_table = coherent_intertwiner_coefficient_table(
        tetrahedra,
        intertwiner_ranges,
        spin_data,
        normal_data,
    )

    first_two_pairs = collect(Iterators.product(
        intertwiner_ranges[1],
        intertwiner_ranges[2],
    ))

    number_of_worker_processes = max(nprocs() - 1, 0)
    number_of_chunks =
        nprocs() == 1 ? 1 :
        min(length(first_two_pairs), number_of_worker_processes * chunks_per_worker)

    chunks = split_into_chunks(first_two_pairs, number_of_chunks)

    if nprocs() == 1
        return su2_bf_intertwiner_chunk_sum(
            simplex,
            spin_data,
            tetrahedra,
            coefficient_table,
            intertwiner_ranges,
            only(chunks),
        )
    end

    partial_sums = pmap(
        chunk -> su2_bf_intertwiner_chunk_sum(
            simplex,
            spin_data,
            tetrahedra,
            coefficient_table,
            intertwiner_ranges,
            chunk,
        ),
        chunks,
    )

    return sum(partial_sums; init = 0.0 + 0.0im)
end

function su2_bf_coherent_lambda_row(lambda::Integer, base_spin, simplex)
    j = lambda * base_spin
    spin_data, normal_data = regular_4simplex_boundary_data(simplex, j)

    closure_error = max_closure_norm(simplex, spin_data, normal_data)

    W = 0.0 + 0.0im
    elapsed = @elapsed W = su2_bf_coherent_amplitude_parallel_over_intertwiners(
        simplex,
        spin_data,
        normal_data,
    )

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
        warmup_su2_bf_worker(base_spin, simplex)
    else
        @sync for worker in workers()
            @async remotecall_fetch(
                warmup_su2_bf_worker,
                worker,
                base_spin,
                simplex,
            )
        end
    end

    return nothing
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

function checkpoint_su2_bf_scan(scan, csv_file, plot_file)
    sorted_scan = sort(scan; by = row -> row.lambda)
    save_su2_bf_scan_csv(csv_file, sorted_scan)
    plot_su2_bf_scan(sorted_scan, plot_file)

    return sorted_scan
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

function print_scan_row(row)
    println(
        row.lambda, "       ",
        Float64(row.j), "     ",
        abs(row.W), "    ",
        angle(row.W), "    ",
        real(row.scaled_W), "    ",
        imag(row.scaled_W), "    ",
        row.elapsed,
    )
    flush(stdout)
end

function run_su2_bf_scan_streaming(lambda_values, csv_file, plot_file)
    scan = []

    println()
    println("lambda  j       |W|                 arg(W)       Re(lambda^6 W)       Im(lambda^6 W)       time")

    for lambda in lambda_values
        row = su2_bf_coherent_lambda_row(lambda, base_spin, simplex)
        push!(scan, row)
        checkpoint_su2_bf_scan(scan, csv_file, plot_file)
        print_scan_row(row)
    end

    return checkpoint_su2_bf_scan(scan, csv_file, plot_file)
end

println("SU(2) BF vertex amplitude with coherent boundary state")
println("vertex backend = fast 6j local vertex + fast 6j recoupling")
println("simplex        = ", simplex)
println("base spin      = ", base_spin)
println("lambda values  = ", lambda_values)
println("processes      = ", nprocs(), " total, ", max(nprocs() - 1, 0), " worker process(es)")
println("parallelism    = over intertwiner chunks inside each lambda")
println("worker launch  = ", get(ENV, "SU2_BF_NWORKERS", "0") == "0" ? "Julia command line / existing workers" : "script addprocs(...; topology=:master_worker)")
println("chunks/worker  = ", chunks_per_worker)
println("cutoff         = none; this is a fixed-boundary coherent vertex amplitude")
println()

println("Warming up workers...")
warmup_su2_bf_workers()

output_dir = joinpath(@__DIR__, "outputs")
mkpath(output_dir)

csv_file = joinpath(output_dir, "su2_bf_coherent_vertex_lambda_scan.csv")
plot_file = joinpath(output_dir, "su2_bf_coherent_lambda6_scan.png")

println("Running lambda scan...")
println("Each completed lambda is checkpointed immediately.")
println("checkpoint data = ", csv_file)
println("checkpoint plot = ", plot_file)

scan = run_su2_bf_scan_streaming(lambda_values, csv_file, plot_file)

println()
println("Expected SU(2) BF large-spin structure:")
println("  W(lambda) ~ lambda^(-6) cos(lambda S_BF + constant phase)")
println("To see the cosine sign oscillation, use Re(lambda^6 W) or Im(lambda^6 W).")
println("|lambda^6 W| is useful as a magnitude/envelope check, but it removes signs.")
println()
println("saved data     = ", csv_file)
println("saved plot     = ", plot_file)
