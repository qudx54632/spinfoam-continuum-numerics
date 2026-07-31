using LinearAlgebra
using Statistics
using Printf
using Plots

# Post-processing check only.
#
# The CSV data are assumed to be produced with the conventions used by the
# code: normalized four-valent intertwiners, the matching coherent-boundary
# coefficients, and the matching 15j vertex amplitude.  This script does not
# convert the coherent coefficient or the 15j convention.
#
# The only idea borrowed from the asymptotic comparison in arXiv:1903.12624 is
# to remove a fitted global coherent-state phase before looking for the Regge
# cosine:
#
#     lambda^6 W_code(lambda)
#       ≈ exp(i (lambda Phi_c + theta0))
#         A cos(5 lambda acos(-1/4) + phi).
#
# Phi_c, theta0, A and phi are therefore convention-dependent fit parameters.

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
    "SU2_BF_GLOBAL_PHASE_PLOT",
    joinpath(csv_dir, csv_stem * "_global_phase_fit.png"),
)

function read_lambda6_scan(filename)
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
                Z = parse(Float64, data["real_lambda6_W"]) +
                    im * parse(Float64, data["imag_lambda6_W"]),
            ),
        )
    end

    rows = filter(row -> row.lambda >= lambda_min, rows)
    sort!(rows; by = row -> row.lambda)

    length(rows) >= 3 ||
        error("need at least three data points to fit the global phase")

    return rows
end

function fit_cosine_with_global_phase(lambdas, values, phi_c, theta0)
    omega = 5 * acos(-1 / 4)
    dephased = exp.(-im .* (phi_c .* lambdas .+ theta0)) .* values
    design_matrix = hcat(cos.(omega .* lambdas), sin.(omega .* lambdas))
    coefficients = design_matrix \ real.(dephased)

    fitted_real = design_matrix * coefficients
    fitted_values = exp.(im .* (phi_c .* lambdas .+ theta0)) .* fitted_real

    residuals = values .- fitted_values
    rmse = sqrt(mean(abs2.(residuals)))
    imaginary_rmse = sqrt(mean(imag.(dephased) .^ 2))

    cosine_coefficient, sine_coefficient = coefficients
    amplitude = hypot(cosine_coefficient, sine_coefficient)
    cosine_offset = atan(-sine_coefficient, cosine_coefficient)

    return (
        phi_c = phi_c,
        theta0 = theta0,
        amplitude = amplitude,
        cosine_offset = cosine_offset,
        dephased = dephased,
        fitted_real = fitted_real,
        fitted_values = fitted_values,
        rmse = rmse,
        imaginary_rmse = imaginary_rmse,
    )
end

function search_global_phase(lambdas, values)
    function search(phi_center, phi_width, theta_center, theta_width, n_phi, n_theta)
        best = fit_cosine_with_global_phase(
            lambdas,
            values,
            phi_center - phi_width,
            theta_center - theta_width,
        )

        for phi_c in range(phi_center - phi_width, phi_center + phi_width; length = n_phi),
            theta0 in range(theta_center - theta_width, theta_center + theta_width; length = n_theta)

            candidate = fit_cosine_with_global_phase(lambdas, values, phi_c, theta0)

            if candidate.rmse < best.rmse
                best = candidate
            end
        end

        return best
    end

    best = search(0.0, pi, 0.0, pi, 721, 721)
    phi_width = 2pi / 720
    theta_width = 2pi / 720

    for _ in 1:5
        best = search(best.phi_c, 4phi_width, best.theta0, 4theta_width, 401, 401)
        phi_width /= 50
        theta_width /= 50
    end

    return best
end

rows = read_lambda6_scan(csv_file)
lambdas = [row.lambda for row in rows]
values = [row.Z for row in rows]

fit = search_global_phase(lambdas, values)

complex_variance = sum(abs2.(values .- mean(values)))
complex_r_squared =
    complex_variance == 0 ? NaN :
    1 - sum(abs2.(values .- fit.fitted_values)) / complex_variance

real_dephased = real.(fit.dephased)
real_variance = sum((real_dephased .- mean(real_dephased)) .^ 2)
real_r_squared =
    real_variance == 0 ? NaN :
    1 - sum((real_dephased .- fit.fitted_real) .^ 2) / real_variance

println("SU(2) BF coherent vertex global-phase fit")
println("amplitude convention = code convention, no coefficient conversion")
println("CSV file          = ", csv_file)
println("lambda range used = ", first(lambdas), ":", last(lambdas))
println()
@printf("Phi_c             = %.15f\n", fit.phi_c)
@printf("theta0            = %.15f\n", fit.theta0)
@printf("A                 = %.15g\n", fit.amplitude)
@printf("cosine offset phi = %.15f\n", fit.cosine_offset)
@printf("complex RMSE      = %.6g\n", fit.rmse)
@printf("complex R^2       = %.6f\n", complex_r_squared)
@printf("dephased Im RMSE  = %.6g\n", fit.imaginary_rmse)
@printf("dephased Re R^2   = %.6f\n", real_r_squared)

p1 = plot(
    lambdas,
    real_dephased;
    marker = :circle,
    linewidth = 2,
    label = "Re(dephased lambda^6 W)",
    xlabel = "lambda",
    ylabel = "dephased lambda^6 W",
    title = "SU(2) BF coherent vertex: fitted global phase removed",
)

plot!(
    p1,
    lambdas,
    fit.fitted_real;
    marker = :diamond,
    linestyle = :dash,
    linewidth = 2,
    label = "fitted Regge cosine",
)

p2 = plot(
    lambdas,
    imag.(fit.dephased);
    marker = :circle,
    linewidth = 2,
    label = "Im(dephased lambda^6 W)",
    xlabel = "lambda",
    ylabel = "imaginary part",
)

hline!(p2, [0.0]; color = :black, linestyle = :dot, label = "0")

plot(p1, p2; layout = (2, 1), size = (900, 750))
savefig(plot_file)

println()
println("saved plot        = ", plot_file)
