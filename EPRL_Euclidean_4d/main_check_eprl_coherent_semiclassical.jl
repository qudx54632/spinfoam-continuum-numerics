include(joinpath(@__DIR__, "semiclassical", "semiclassical_vertex_tools.jl"))

using Printf

# Server-safe Euclidean EPRL coherent-boundary single-vertex scan.
#
# This is the EPRL analogue of `main_check_su2_bf_semiclassical.jl`.
# For each fixed lambda, independent Julia child processes split the final
# right-sector intertwiner sum.  The diagonal boundary intertwiners are summed
# locally into coherent-fusion coefficients before the vertex contraction.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

const simplex = (1, 2, 3, 4, 5)

function parse_rational_setting(name, default_value)
    text = get(ENV, name, default_value)

    if occursin("//", text)
        return parse(Rational{Int}, text)
    else
        return parse(Int, text) // 1
    end
end

const gamma = parse_rational_setting("EPRL_GAMMA", "1//3")
const base_spin = parse_rational_setting("EPRL_BASE_SPIN", "3//2")

lambda_start = parse(Int, get(ENV, "EPRL_LAMBDA_START", "1"))
lambda_stop = parse(Int, get(ENV, "EPRL_LAMBDA_STOP", "3"))
lambda_step = parse(Int, get(ENV, "EPRL_LAMBDA_STEP", "1"))
lambda_values = collect(lambda_start:lambda_step:lambda_stop)

# Number of independent Julia child processes used inside each lambda.
# Use this instead of `julia -p N`.
process_workers = parse(Int, get(ENV, "EPRL_NWORKERS", "1"))

const output_dir = joinpath(@__DIR__, "outputs")
const partial_dir = joinpath(output_dir, "eprl_coherent_partial_sums")
const csv_file = joinpath(output_dir, "eprl_coherent_vertex_lambda_scan.csv")
const plot_file = joinpath(output_dir, "eprl_coherent_lambda12_scan.png")

function split_into_chunks(values, number_of_chunks)
    isempty(values) && return []

    number_of_chunks = max(1, min(number_of_chunks, length(values)))
    chunk_size = cld(length(values), number_of_chunks)

    return [
        values[first:min(first + chunk_size - 1, length(values))]
        for first in 1:chunk_size:length(values)
    ]
end

function coherent_boundary_data_for_lambda(lambda)
    j = lambda * base_spin
    spin_data, normal_data = paper_twisted_spike_boundary_data(simplex, j)

    return j, spin_data, normal_data
end

function eprl_boundary_spins_allowed(spin_data, gamma)
    return all(eprl_spin_allowed(j, gamma) for j in values(spin_data))
end

function intertwiner_data_from_tuple(tetrahedra, intertwiners)
    return Dict(
        tetrahedron => i
        for (tetrahedron, i) in zip(tetrahedra, intertwiners)
    )
end

"""
For each boundary tetrahedron tau, precompute

    C_tau(iL,iR) = sum_i c_i(tau) f^i_{iL iR}(tau).

This is the key speedup: the coherent boundary coefficient and the EPRL
fusion map are both local to one tetrahedron, so the diagonal intertwiner
sum can be done before summing the left and right SU(2) BF vertices.
"""
function coherent_fusion_tables(
    tetrahedra,
    spin_data::AbstractDict,
    normal_data::AbstractDict,
    gamma,
)
    diagonal_ranges = [
        allowed_global_intertwiners(tetrahedron, spin_data)
        for tetrahedron in tetrahedra
    ]

    coefficient_table = coherent_intertwiner_coefficient_table(
        tetrahedra,
        diagonal_ranges,
        spin_data,
        normal_data,
    )

    left_spin_data = eprl_left_spin_data(spin_data, gamma)
    right_spin_data = eprl_right_spin_data(spin_data, gamma)

    left_ranges = [
        allowed_global_intertwiners(tetrahedron, left_spin_data)
        for tetrahedron in tetrahedra
    ]
    right_ranges = [
        allowed_global_intertwiners(tetrahedron, right_spin_data)
        for tetrahedron in tetrahedra
    ]

    tables = [Dict() for _ in tetrahedra]

    for (a, tetrahedron) in enumerate(tetrahedra)
        spins = global_intertwiner_spins(tetrahedron, spin_data)

        for i in diagonal_ranges[a]
            c = coefficient_table[(tetrahedron, i)]
            abs(c) == 0.0 && continue

            for (iL, iR) in allowed_eprl_intertwiner_pairs(
                tetrahedron,
                spin_data,
                i,
                gamma,
            )
                f = eprl_fusion_coefficient(spins, i, iL, iR, gamma)
                f == 0.0 && continue

                key = (iL, iR)
                tables[a][key] = get(tables[a], key, 0.0 + 0.0im) + c * f
            end
        end
    end

    return (
        left_spin_data = left_spin_data,
        right_spin_data = right_spin_data,
        left_ranges = left_ranges,
        right_ranges = right_ranges,
        tables = tables,
    )
end

function right_intertwiner_tuple_chunks(lambda, number_of_chunks)
    _, spin_data, _ = coherent_boundary_data_for_lambda(lambda)

    eprl_boundary_spins_allowed(spin_data, gamma) || return []

    tetrahedra = oriented_simplex_tetrahedra(simplex)
    right_spin_data = eprl_right_spin_data(spin_data, gamma)
    right_ranges = [
        allowed_global_intertwiners(tetrahedron, right_spin_data)
        for tetrahedron in tetrahedra
    ]

    any(isempty, right_ranges) && return []

    right_tuples = collect(Iterators.product(right_ranges...))

    return split_into_chunks(right_tuples, number_of_chunks)
end

function coherent_fusion_product(tables, left_tuple, right_tuple)
    product = 1.0 + 0.0im

    for a in 1:length(tables)
        c = get(tables[a], (left_tuple[a], right_tuple[a]), 0.0 + 0.0im)
        abs(c) == 0.0 && return 0.0 + 0.0im
        product *= c
    end

    return product
end

function eprl_rearranged_right_chunk_sum(
    simplex,
    spin_data::AbstractDict,
    normal_data::AbstractDict,
    right_tuple_chunk,
    gamma,
)
    sigma = oriented_four_simplex(simplex)
    tetrahedra = oriented_simplex_tetrahedra(sigma)

    fusion_data = coherent_fusion_tables(
        tetrahedra,
        spin_data,
        normal_data,
        gamma,
    )

    left_vertex_table = simplex_amplitude_table_recoupled_oriented_6j(
        sigma,
        fusion_data.left_spin_data,
        fusion_data.left_ranges,
    )
    right_vertex_table = simplex_amplitude_table_recoupled_oriented_6j(
        sigma,
        fusion_data.right_spin_data,
        fusion_data.right_ranges,
    )

    left_vertices = [
        (left_tuple = left_tuple, vertex = vertex)
        for (left_tuple, vertex) in left_vertex_table
        if vertex != 0.0
    ]

    total = 0.0 + 0.0im

    for right_tuple in right_tuple_chunk
        right_vertex = get(right_vertex_table, right_tuple, 0.0)

        right_vertex == 0.0 && continue

        for left_entry in left_vertices
            coefficient = coherent_fusion_product(
                fusion_data.tables,
                left_entry.left_tuple,
                right_tuple,
            )

            abs(coefficient) == 0.0 && continue
            total += left_entry.vertex * right_vertex * coefficient
        end
    end

    return total
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
    lambda = parse(Int, ENV["EPRL_WORKER_LAMBDA"])
    chunk_id = parse(Int, ENV["EPRL_WORKER_CHUNK_ID"])
    nchunks = parse(Int, ENV["EPRL_WORKER_NCHUNKS"])
    output_file = ENV["EPRL_WORKER_OUTPUT_FILE"]

    _, spin_data, normal_data = coherent_boundary_data_for_lambda(lambda)

    if !eprl_boundary_spins_allowed(spin_data, gamma)
        write_partial_sum(output_file, lambda, chunk_id, nchunks, 0.0 + 0.0im, 0.0)
        return nothing
    end

    chunks = right_intertwiner_tuple_chunks(lambda, nchunks)
    chunk = chunk_id <= length(chunks) ? chunks[chunk_id] : []

    partial_sum = 0.0 + 0.0im
    elapsed = @elapsed partial_sum = eprl_rearranged_right_chunk_sum(
        simplex,
        spin_data,
        normal_data,
        chunk,
        gamma,
    )

    write_partial_sum(output_file, lambda, chunk_id, nchunks, partial_sum, elapsed)

    return nothing
end

if get(ENV, "EPRL_PARTIAL_WORKER", "0") == "1"
    run_partial_sum_worker()
    exit()
end

using Plots

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

function save_eprl_scan_csv(filename, scan)
    open(filename, "w") do io
        println(
            io,
            "lambda,j,allowed,abs_W,arg_W,real_W,imag_W,real_lambda12_W,imag_lambda12_W,abs_lambda12_W,S_int,S_ext,closure_error,elapsed_seconds",
        )

        for row in scan
            println(
                io,
                row.lambda, ",",
                Float64(row.j), ",",
                row.allowed, ",",
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

function plot_eprl_scan(scan, plot_file)
    sorted_scan = sort(scan; by = row -> row.lambda)
    x = [row.lambda for row in sorted_scan]
    phase = dephasing_phase(sorted_scan)

    dephased_values = [
        exp(-im * phase * row.lambda) * row.scaled_W
        for row in sorted_scan
    ]
    numerical_real = real.(dephased_values)
    numerical_imag = imag.(dephased_values)

    p_real = scatter(
        x,
        numerical_real;
        color = :darkred,
        markerstrokecolor = :darkred,
        marker = :circle,
        label = "Numerical data",
        xlabel = "lambda",
        ylabel = "Re[dephased lambda^12 W]",
        title = "Euclidean EPRL coherent vertex: real part",
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
        ylabel = "Im[dephased lambda^12 W]",
        title = "Imaginary part",
        yticks = symmetric_ticks(numerical_imag),
        yformatter = y -> @sprintf("%.1e", y),
    )
    hline!(p_imag, [0.0]; color = :black, linewidth = 1, label = "")

    plot(p_real, p_imag; layout = (1, 2), size = (1200, 500))
    savefig(plot_file)
end

function checkpoint_eprl_scan(scan)
    sorted_scan = sort(scan; by = row -> row.lambda)
    save_eprl_scan_csv(csv_file, sorted_scan)
    plot_eprl_scan(sorted_scan, plot_file)

    return sorted_scan
end

function worker_command(lambda, chunk_id, nchunks, output_file)
    command = `$(Base.julia_cmd()) --startup-file=no $(abspath(PROGRAM_FILE))`

    return addenv(
        command,
        "EPRL_PARTIAL_WORKER" => "1",
        "EPRL_WORKER_LAMBDA" => string(lambda),
        "EPRL_WORKER_CHUNK_ID" => string(chunk_id),
        "EPRL_WORKER_NCHUNKS" => string(nchunks),
        "EPRL_WORKER_OUTPUT_FILE" => output_file,
    )
end

function run_lambda_with_process_chunks(lambda)
    mkpath(partial_dir)

    j, spin_data, normal_data = coherent_boundary_data_for_lambda(lambda)
    closure_error = max_closure_norm(simplex, spin_data, normal_data)
    allowed = eprl_boundary_spins_allowed(spin_data, gamma)

    S_int = regular_4simplex_regge_action(j, gamma; dihedral = :interior)
    S_ext = regular_4simplex_regge_action(j, gamma; dihedral = :exterior)

    if !allowed
        W = 0.0 + 0.0im
        return (
            lambda = lambda,
            j = j,
            allowed = allowed,
            W = W,
            scaled_W = Float64(lambda)^12 * W,
            S_int = S_int,
            S_ext = S_ext,
            closure_error = closure_error,
            elapsed = 0.0,
        )
    end

    chunks = right_intertwiner_tuple_chunks(lambda, process_workers)
    nchunks = length(chunks)

    nchunks == 0 && begin
        W = 0.0 + 0.0im
        return (
            lambda = lambda,
            j = j,
            allowed = allowed,
            W = W,
            scaled_W = Float64(lambda)^12 * W,
            S_int = S_int,
            S_ext = S_ext,
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

    return (
        lambda = lambda,
        j = j,
        allowed = allowed,
        W = W,
        scaled_W = Float64(lambda)^12 * W,
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
        row.allowed, "     ",
        abs(row.W), "    ",
        angle(row.W), "    ",
        real(row.scaled_W), "    ",
        imag(row.scaled_W), "    ",
        row.elapsed,
    )
    flush(stdout)
end

function print_eprl_check_summary(scan)
    isempty(scan) && return nothing

    sorted_scan = sort(scan; by = row -> row.lambda)
    phase = dephasing_phase(sorted_scan)

    println()
    println("Paper-style coherent EPRL check:")
    println("  plotted data = Re/Im[exp(-i Phi_c lambda) lambda^12 W(lambda)]")
    println("  Phi_c        = ", phase)
    println("  S_int        = gamma * 10 * j * acos(1/4)")
    println("  S_ext        = gamma * 10 * j * acos(-1/4)")

    return nothing
end

println("Euclidean EPRL vertex amplitude with coherent boundary state")
println("vertex backend = EPRL fusion map + fast 6j SU(2) BF vertices")
println("summation      = local coherent-fusion tables C_tau(iL,iR)")
println("parallelism    = independent Julia processes over right-sector chunks")
println("simplex        = ", simplex)
println("gamma          = ", gamma)
println("base spin      = ", base_spin)
println("boundary data  = exact twisted-spike normals")
println("lambda values  = ", lambda_values)
println("chunk workers  = ", process_workers, " per lambda")
println("cutoff         = none; this is a fixed-boundary coherent vertex amplitude")
println("plot convention = dephase lambda^12 W and keep signs")
println()

mkpath(output_dir)
mkpath(partial_dir)

println("Each completed lambda is checkpointed immediately.")
println("checkpoint data = ", csv_file)
println("checkpoint plot = ", plot_file)
println()
println("lambda  j       allowed  |W|                 arg(W)       Re(lambda^12 W)      Im(lambda^12 W)      time")

scan = []

for lambda in lambda_values
    row = run_lambda_with_process_chunks(lambda)
    push!(scan, row)
    checkpoint_eprl_scan(scan)
    print_scan_row(row)
end

println()
println("Expected Euclidean EPRL large-spin structure:")
println("  W(lambda) ~ lambda^(-12) times Regge-frequency oscillatory terms")
println("The checkpoint plot removes the coherent-state global phase and keeps signs.")
print_eprl_check_summary(scan)
println()
println("saved data     = ", csv_file)
println("saved plot     = ", plot_file)
