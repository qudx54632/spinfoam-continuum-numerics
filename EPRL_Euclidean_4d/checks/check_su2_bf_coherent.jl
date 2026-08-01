using Printf

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
using Plots

default_csv = joinpath(
    dirname(@__DIR__),
    "outputs",
    "su2_bf_coherent_vertex_lambda_scan.csv",
)

csv_file = get(ENV, "SU2_BF_SCAN_CSV", default_csv)
csv_dir = dirname(csv_file)
csv_stem = splitext(basename(csv_file))[1]

plot_file = get(
    ENV,
    "SU2_BF_COHERENT_CHECK_PLOT",
    get(
        ENV,
        "SU2_BF_1708_STYLE_PLOT",
        joinpath(csv_dir, csv_stem * "_coherent_check.png"),
    ),
)

function read_scan_csv(filename)
    lines = readlines(filename)
    isempty(lines) && error("empty CSV file: $filename")

    header = split(lines[1], ",")
    rows = NamedTuple[]

    for line in lines[2:end]
        isempty(strip(line)) && continue

        fields = split(line, ",")
        data = Dict(header .=> fields)

        push!(
            rows,
            (
                lambda = parse(Float64, data["lambda"]),
                lambda6_W = parse(Float64, data["real_lambda6_W"]) +
                    im * parse(Float64, data["imag_lambda6_W"]),
            ),
        )
    end

    # The scan parameter is lambda.  The `j` column is only the
    # corresponding boundary spin (here j=lambda/2); it is not the
    # horizontal variable and is not used for the lambda^6 rescaling.
    sort!(rows; by = row -> row.lambda)
    return rows
end

rows = read_scan_csv(csv_file)
x = [row.lambda for row in rows]

# The current coherent-state convention carries a linear global phase.
# Removing it is the analogue of the e^{-i lambda S_R} step in Fig. 4.
function dephasing_phase(rows)
    for row in rows
        row.lambda == 0 && continue
        abs(row.lambda6_W) == 0 && continue

        return angle(row.lambda6_W) / row.lambda
    end

    return 0.0
end

phase_per_lambda = dephasing_phase(rows)
numerical_data = [
    exp(-im * phase_per_lambda * row.lambda) * row.lambda6_W
    for row in rows
]
numerical_real = real.(numerical_data)
numerical_imag = imag.(numerical_data)

# Eq. (LO2) in arXiv:1708.01727, adapted to the current base spin
# j_ab = lambda/2.  The paper writes the equilateral Regge frequency with
# the exterior angle acos(-1/4).  The present coherent-state convention is
# naturally compared with the interior angle acos(1/4), after the numerical
# global phase above has been removed.
analytic_dihedral = get(ENV, "SU2_BF_ANALYTIC_DIHEDRAL", "interior")
analytic_dihedral in ("interior", "exterior") ||
    error("SU2_BF_ANALYTIC_DIHEDRAL must be either interior or exterior")

theta = analytic_dihedral == "interior" ? acos(1 / 4) : acos(-1 / 4)
omega_lambda = 5 * theta

det_unit =
    (1618200 - im * 316712 * sqrt(15)) / 177147
det_base = (1 / 2)^12 * det_unit
paper_amplitude =
    (2pi)^6 * 2^4 / (4pi)^8 / sqrt(abs(det_base))
paper_phase =
    analytic_dihedral == "interior" ? angle(det_base) / 2 : -angle(det_base) / 2
paper_points =
    paper_amplitude .* cos.(omega_lambda .* x .+ paper_phase)

function symmetric_ticks(values)
    scale = maximum(abs.(values))
    scale = scale == 0 ? 1.0 : scale
    step = 10.0 ^ floor(log10(scale))
    top = ceil(scale / step) * step
    return collect(range(-top, top; length = 5))
end

p_real = scatter(
    x,
    numerical_real;
    color = :darkred,
    markerstrokecolor = :darkred,
    markersize = 5,
    label = "Numerical Data",
    xlabel = "lambda",
    ylabel = "Re[dephased lambda^6 A_v]",
    title = "Real part, $(analytic_dihedral) angle convention",
)

scatter!(
    p_real,
    x,
    paper_points;
    color = :white,
    markerstrokecolor = :deepskyblue3,
    markersize = 5,
    label = "Analytic samples",
)

hline!(p_real, [0.0]; color = :black, linewidth = 1, label = "")

p_imag = scatter(
    x,
    numerical_imag;
    color = :darkblue,
    markerstrokecolor = :darkblue,
    markersize = 5,
    label = "Numerical Data",
    xlabel = "lambda",
    ylabel = "Im[dephased lambda^6 A_v]",
    title = "Imaginary part",
    yticks = symmetric_ticks(numerical_imag),
    yformatter = y -> @sprintf("%.1e", y),
)

hline!(p_imag, [0.0]; color = :black, linewidth = 1, label = "")

plot(p_real, p_imag; layout = (1, 2), size = (1200, 500))
savefig(plot_file)

println("CSV file      = ", csv_file)
println("plot          = ", plot_file)
println()
@printf("phase_per_lambda     = %.15f\n", phase_per_lambda)
println("analytic dihedral    = ", analytic_dihedral)
@printf("analytic omega       = %.15f\n", omega_lambda)
@printf("paper LO2 A   = %.15g\n", paper_amplitude)
@printf("paper LO2 phi = %.15f\n", paper_phase)
