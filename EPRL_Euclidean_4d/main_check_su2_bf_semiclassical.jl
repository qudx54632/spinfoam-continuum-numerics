include(joinpath(@__DIR__, "semiclassical", "semiclassical_vertex_tools.jl"))

using Printf

# Server-safe SU(2) BF coherent-boundary semiclassical scan.
#
# This script does not use Julia Distributed.  Instead, for each fixed lambda
# it launches independent Julia child processes, each evaluating one chunk of
# the boundary-intertwiner sum.  This avoids the TCP worker-credential errors
# that can appear with `julia -p N` on some servers.

const simplex = (1, 2, 3, 4, 5)
const base_spin = 1 // 2

lambda_start = parse(Int, get(ENV, "SU2_BF_LAMBDA_START", "1"))
lambda_stop = parse(Int, get(ENV, "SU2_BF_LAMBDA_STOP", "70"))
lambda_values = collect(lambda_start:lambda_stop)

# Number of independent Julia child processes used inside each lambda.
# Use this instead of `julia -p N`.
process_workers = parse(Int, get(ENV, "SU2_BF_NWORKERS", "1"))

const output_dir = joinpath(@__DIR__, "outputs")
const partial_dir = joinpath(output_dir, "su2_bf_partial_sums")
const csv_file = joinpath(output_dir, "su2_bf_coherent_vertex_lambda_scan.csv")
const plot_file = joinpath(output_dir, "su2_bf_coherent_lambda6_scan.png")

function split_into_chunks(values, number_of_chunks)
    isempty(values) && return []

    number_of_chunks = max(1, min(number_of_chunks, length(values)))
    chunk_size = cld(length(values), number_of_chunks)

    return [
        values[first:min(first + chunk_size - 1, length(values))]
        for first in 1:chunk_size:length(values)
    ]
end

function su2_bf_intertwiner_chunk_sum(
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

function coherent_boundary_data_for_lambda(lambda)
    j = lambda * base_spin
    spin_data, normal_data = paper_twisted_spike_boundary_data(simplex, j)

    return j, spin_data, normal_data
end

function first_two_intertwiner_pair_chunks(lambda, number_of_chunks)
    _, spin_data, _ = coherent_boundary_data_for_lambda(lambda)
    tetrahedra = oriented_simplex_tetrahedra(simplex)
    intertwiner_ranges = [
        allowed_global_intertwiners(tetrahedron, spin_data)
        for tetrahedron in tetrahedra
    ]

    any(isempty, intertwiner_ranges) && return []

    first_two_pairs = collect(Iterators.product(
        intertwiner_ranges[1],
        intertwiner_ranges[2],
    ))

    return split_into_chunks(first_two_pairs, number_of_chunks)
end

function partial_sum_filename(lambda, chunk_id)
    return joinpath(
        partial_dir,
        "lambda_$(lambda)_chunk_$(chunk_id).csv",
    )
end

function write_partial_sum(filename, lambda, chunk_id, nchunks, partial_sum, elapsed)
    open(filename, "w") do io
        println(io, "lambda,chunk_id,nchunks,real_partial,imag_partial,elapsed_seconds")
        println(
            io,
            lambda, ",",
            chunk_id, ",",
            nchunks, ",",
            real(partial_sum), ",",
            imag(partial_sum), ",",
            elapsed,
        )
    end
end

function read_partial_sum(filename)
    lines = readlines(filename)
    length(lines) >= 2 ||
        throw(ArgumentError("partial sum file is empty: $filename"))

    fields = split(lines[2], ",")
    length(fields) == 6 ||
        throw(ArgumentError("partial sum file has wrong format: $filename"))

    return (
        lambda = parse(Int, fields[1]),
        chunk_id = parse(Int, fields[2]),
        nchunks = parse(Int, fields[3]),
        partial_sum = parse(Float64, fields[4]) + im * parse(Float64, fields[5]),
        elapsed = parse(Float64, fields[6]),
    )
end

function run_partial_sum_worker()
    lambda = parse(Int, ENV["SU2_BF_WORKER_LAMBDA"])
    chunk_id = parse(Int, ENV["SU2_BF_WORKER_CHUNK_ID"])
    nchunks = parse(Int, ENV["SU2_BF_WORKER_NCHUNKS"])
    output_file = ENV["SU2_BF_WORKER_OUTPUT_FILE"]

    _, spin_data, normal_data = coherent_boundary_data_for_lambda(lambda)
    tetrahedra = oriented_simplex_tetrahedra(simplex)

    intertwiner_ranges = [
        allowed_global_intertwiners(tetrahedron, spin_data)
        for tetrahedron in tetrahedra
    ]

    any(isempty, intertwiner_ranges) && begin
        write_partial_sum(output_file, lambda, chunk_id, nchunks, 0.0 + 0.0im, 0.0)
        return nothing
    end

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

    chunks = split_into_chunks(first_two_pairs, nchunks)
    chunk =
        chunk_id <= length(chunks) ? chunks[chunk_id] :
        Tuple{eltype(first_two_pairs)}[]

    partial_sum = 0.0 + 0.0im
    elapsed = @elapsed partial_sum = su2_bf_intertwiner_chunk_sum(
        simplex,
        spin_data,
        tetrahedra,
        coefficient_table,
        intertwiner_ranges,
        chunk,
    )

    write_partial_sum(output_file, lambda, chunk_id, nchunks, partial_sum, elapsed)

    return nothing
end

if get(ENV, "SU2_BF_PARTIAL_WORKER", "0") == "1"
    run_partial_sum_worker()
    exit()
end

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots

function su2_bf_analytic_parameters()
    dihedral = get(ENV, "SU2_BF_ANALYTIC_DIHEDRAL", "interior")
    dihedral in ("interior", "exterior") ||
        error("SU2_BF_ANALYTIC_DIHEDRAL must be either interior or exterior")

    theta = dihedral == "interior" ? acos(1 / 4) : acos(-1 / 4)
    omega = 5 * theta

    det_unit = (1618200 - im * 316712 * sqrt(15)) / 177147
    det_base = (1 / 2)^12 * det_unit
    amplitude = (2pi)^6 * 2^4 / (4pi)^8 / sqrt(abs(det_base))
    phase = dihedral == "interior" ? angle(det_base) / 2 : -angle(det_base) / 2

    return (
        dihedral = dihedral,
        omega = omega,
        amplitude = amplitude,
        phase = phase,
    )
end

function dephasing_phase(sorted_scan)
    for row in sorted_scan
        row.lambda == 0 && continue
        abs(row.scaled_W) == 0 && continue

        return angle(row.scaled_W) / row.lambda
    end

    return 0.0
end

function symmetric_ticks(values)
    isempty(values) && return [-1.0, 0.0, 1.0]

    scale = maximum(abs.(values))
    scale = scale == 0 ? 1.0 : scale
    step = 10.0 ^ floor(log10(scale))
    top = ceil(scale / step) * step

    return collect(range(-top, top; length = 5))
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
    sorted_scan = sort(scan; by = row -> row.lambda)
    x = [row.lambda for row in sorted_scan]
    phase = dephasing_phase(sorted_scan)

    dephased_values = [
        exp(-im * phase * row.lambda) * row.scaled_W
        for row in sorted_scan
    ]
    numerical_real = real.(dephased_values)
    numerical_imag = imag.(dephased_values)

    analytic = su2_bf_analytic_parameters()
    analytic_points =
        analytic.amplitude .* cos.(analytic.omega .* x .+ analytic.phase)

    p_real = scatter(
        x,
        numerical_real;
        color = :darkred,
        markerstrokecolor = :darkred,
        marker = :circle,
        label = "Numerical data",
        xlabel = "lambda",
        ylabel = "Re[dephased lambda^6 W]",
        title = "Real part, $(analytic.dihedral) angle convention",
    )

    scatter!(
        p_real,
        x,
        analytic_points;
        color = :white,
        markerstrokecolor = :deepskyblue3,
        marker = :circle,
        label = "Analytic samples",
    )

    hline!(p_real, [0.0]; color = :black, linewidth = 1, label = "")

    p_imag = scatter(
        x,
        numerical_imag;
        color = :darkblue,
        markerstrokecolor = :darkblue,
        marker = :circle,
        label = "Numerical data",
        xlabel = "lambda",
        ylabel = "Im[dephased lambda^6 W]",
        title = "Imaginary part",
        yticks = symmetric_ticks(numerical_imag),
        yformatter = y -> @sprintf("%.1e", y),
    )

    hline!(p_imag, [0.0]; color = :black, linewidth = 1, label = "")

    plot(p_real, p_imag; layout = (1, 2), size = (1200, 500))
    savefig(plot_file)
end

function print_su2_bf_check_summary(scan)
    sorted_scan = sort(scan; by = row -> row.lambda)
    phase = dephasing_phase(sorted_scan)
    analytic = su2_bf_analytic_parameters()

    println()
    println("Paper-style coherent SU(2) BF check:")
    println("  plotted data       = Re/Im[exp(-i Phi_c lambda) lambda^6 W(lambda)]")
    println("  Phi_c              = ", phase)
    println("  analytic dihedral  = ", analytic.dihedral)
    println("  analytic samples   = A cos(omega lambda + phi)")
    println("  A                  = ", analytic.amplitude)
    println("  omega              = ", analytic.omega)
    println("  phi                = ", analytic.phase)

    return nothing
end

function checkpoint_su2_bf_scan(scan)
    sorted_scan = sort(scan; by = row -> row.lambda)
    save_su2_bf_scan_csv(csv_file, sorted_scan)
    plot_su2_bf_scan(sorted_scan, plot_file)

    return sorted_scan
end

function worker_command(lambda, chunk_id, nchunks, output_file)
    command = `$(Base.julia_cmd()) --startup-file=no $(abspath(PROGRAM_FILE))`

    return addenv(
        command,
        "SU2_BF_PARTIAL_WORKER" => "1",
        "SU2_BF_WORKER_LAMBDA" => string(lambda),
        "SU2_BF_WORKER_CHUNK_ID" => string(chunk_id),
        "SU2_BF_WORKER_NCHUNKS" => string(nchunks),
        "SU2_BF_WORKER_OUTPUT_FILE" => output_file,
    )
end

function run_lambda_with_process_chunks(lambda)
    mkpath(partial_dir)

    _, spin_data, normal_data = coherent_boundary_data_for_lambda(lambda)
    closure_error = max_closure_norm(simplex, spin_data, normal_data)

    chunks = first_two_intertwiner_pair_chunks(lambda, process_workers)
    nchunks = length(chunks)

    nchunks == 0 && begin
        j = lambda * base_spin
        W = 0.0 + 0.0im
        return (
            lambda = lambda,
            j = j,
            W = W,
            scaled_W = lambda^6 * W,
            S_int = regular_4simplex_su2_bf_regge_action(j; dihedral = :interior),
            S_ext = regular_4simplex_su2_bf_regge_action(j; dihedral = :exterior),
            closure_error = closure_error,
            elapsed = 0.0,
        )
    end

    output_files = [partial_sum_filename(lambda, chunk_id) for chunk_id in 1:nchunks]

    for output_file in output_files
        isfile(output_file) && rm(output_file)
    end

    elapsed_total = @elapsed begin
        processes = [
            run(
                worker_command(lambda, chunk_id, nchunks, output_files[chunk_id]);
                wait = false,
            )
            for chunk_id in 1:nchunks
        ]

        for process in processes
            wait(process)
            success(process) || error("a partial worker failed for lambda=$lambda")
        end
    end

    partial_rows = [read_partial_sum(output_file) for output_file in output_files]
    W = sum(row.partial_sum for row in partial_rows; init = 0.0 + 0.0im)

    j = lambda * base_spin
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
        elapsed = elapsed_total,
    )
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

println("SU(2) BF vertex amplitude with coherent boundary state")
println("vertex backend = fast 6j local vertex + fast 6j recoupling")
println("parallelism    = independent Julia processes over intertwiner chunks")
println("simplex        = ", simplex)
println("base spin      = ", base_spin)
println("boundary data  = paper twisted-spike normals, arXiv:1903.12624 Table 2")
println("lambda values  = ", lambda_values)
println("chunk workers  = ", process_workers, " per lambda")
println("cutoff         = none; this is a fixed-boundary coherent vertex amplitude")
println("plot convention = dephase lambda^6 W and compare with analytic sample points")
println()

mkpath(output_dir)
mkpath(partial_dir)

println("Each completed lambda is checkpointed immediately.")
println("checkpoint data = ", csv_file)
println("checkpoint plot = ", plot_file)
println()
println("lambda  j       |W|                 arg(W)       Re(lambda^6 W)       Im(lambda^6 W)       time")

scan = []

for lambda in lambda_values
    row = run_lambda_with_process_chunks(lambda)
    push!(scan, row)
    checkpoint_su2_bf_scan(scan)
    print_scan_row(row)
end

println()
println("Expected SU(2) BF large-spin structure:")
println("  W(lambda) ~ lambda^(-6) cos(lambda S_BF + constant phase)")
println("The checkpoint plot removes the coherent-state global phase and keeps signs.")
print_su2_bf_check_summary(scan)
println()
println("saved data     = ", csv_file)
println("saved plot     = ", plot_file)
