using LinearAlgebra
using Statistics
using Printf
using Plots

default_csv = joinpath(
    dirname(@__DIR__),
    "outputs",
    "su2_bf_coherent_vertex_lambda_scan.csv",
)

csv_file = get(ENV, "SU2_BF_SCAN_CSV", default_csv)
lambda_min = parse(Float64, get(ENV, "SU2_BF_FIT_LAMBDA_MIN", "-Inf"))

csv_dir = dirname(csv_file)
csv_stem = splitext(basename(csv_file))[1]
plot_file = get(
    ENV,
    "SU2_BF_FIT_PLOT",
    joinpath(csv_dir, csv_stem * "_regge_fit.png"),
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
                real_lambda6_W = parse(Float64, data["real_lambda6_W"]),
                imag_lambda6_W = parse(Float64, data["imag_lambda6_W"]),
                abs_lambda6_W = parse(Float64, data["abs_lambda6_W"]),
                S_ext = parse(Float64, data["S_ext"]),
                S_int = parse(Float64, data["S_int"]),
            ),
        )
    end

    return rows
end

function fit_cosine(phases, values; offset = false)
    design_matrix =
        offset ? hcat(cos.(phases), sin.(phases), ones(length(phases))) :
        hcat(cos.(phases), sin.(phases))

    coefficients = design_matrix \ values
    fitted_values = design_matrix * coefficients

    cosine_coefficient, sine_coefficient = coefficients[1], coefficients[2]
    amplitude = hypot(cosine_coefficient, sine_coefficient)
    phase = atan(-sine_coefficient, cosine_coefficient)
    constant_offset = offset ? coefficients[3] : 0.0

    residuals = values .- fitted_values
    rmse = sqrt(mean(residuals .^ 2))
    total_variance = sum((values .- mean(values)) .^ 2)
    r_squared =
        total_variance == 0 ? NaN :
        1 - sum(residuals .^ 2) / total_variance

    return (
        amplitude = amplitude,
        phase = phase,
        offset = constant_offset,
        fitted_values = fitted_values,
        rmse = rmse,
        r_squared = r_squared,
    )
end

scan = filter(row -> row.lambda >= lambda_min, read_scan_csv(csv_file))
sort!(scan; by = row -> row.lambda)

length(scan) >= 2 ||
    error("need at least two data points for the Regge cosine fit")

lambdas = [row.lambda for row in scan]
real_scaled = [row.real_lambda6_W for row in scan]
imag_scaled = [row.imag_lambda6_W for row in scan]
abs_scaled = [row.abs_lambda6_W for row in scan]
S_ext = [row.S_ext for row in scan]
S_int = [row.S_int for row in scan]

real_fit_ext = fit_cosine(S_ext, real_scaled)
imag_fit_ext = fit_cosine(S_ext, imag_scaled)
real_fit_int = fit_cosine(S_int, real_scaled)
imag_fit_int = fit_cosine(S_int, imag_scaled)

println("Regge cosine fit for SU(2) BF coherent vertex data")
println("CSV file = ", csv_file)
println("lambda range used = ", first(lambdas), ":", last(lambdas))
println()

@printf(
    "Re(lambda^6 W) ≈ %.12g cos(S_ext + %.12g),  R² = %.6f, RMSE = %.4g\n",
    real_fit_ext.amplitude,
    real_fit_ext.phase,
    real_fit_ext.r_squared,
    real_fit_ext.rmse,
)

@printf(
    "Im(lambda^6 W) ≈ %.12g cos(S_ext + %.12g),  R² = %.6f, RMSE = %.4g\n",
    imag_fit_ext.amplitude,
    imag_fit_ext.phase,
    imag_fit_ext.r_squared,
    imag_fit_ext.rmse,
)

println()
println("For comparison, using the interior-angle frequency gives:")

@printf(
    "Re(lambda^6 W) with S_int: A = %.12g, phi = %.12g, R² = %.6f\n",
    real_fit_int.amplitude,
    real_fit_int.phase,
    real_fit_int.r_squared,
)

@printf(
    "Im(lambda^6 W) with S_int: A = %.12g, phi = %.12g, R² = %.6f\n",
    imag_fit_int.amplitude,
    imag_fit_int.phase,
    imag_fit_int.r_squared,
)

p1 = plot(
    lambdas,
    real_scaled;
    marker = :circle,
    label = "Re(lambda^6 W)",
    ylabel = "lambda^6 W",
)

plot!(
    p1,
    lambdas,
    real_fit_ext.fitted_values;
    linestyle = :dash,
    label = "Re fit",
)

plot!(
    p1,
    lambdas,
    imag_scaled;
    marker = :diamond,
    label = "Im(lambda^6 W)",
)

plot!(
    p1,
    lambdas,
    imag_fit_ext.fitted_values;
    linestyle = :dashdot,
    label = "Im fit",
)

p2 = plot(
    lambdas,
    abs_scaled;
    marker = :circle,
    xlabel = "lambda",
    ylabel = "|lambda^6 W|",
    label = "|lambda^6 W|",
)

plot(p1, p2; layout = (2, 1), size = (900, 700))
savefig(plot_file)

println()
println("saved fit plot = ", plot_file)
