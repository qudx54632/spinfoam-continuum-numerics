using Random

include(joinpath(@__DIR__, "..", "load_eprl.jl"))

const validation_tolerance = 1.0e-10

function local_intertwiner_ranges(J)
    return [
        [
            i
            for i in half_integer_spins(min(spins[1] + spins[2], spins[3] + spins[4]))
            if admissible_4valent(spins, i)
        ]
        for spins in intertwiner_spin_lists(J)
    ]
end

function local_intertwiner_tuples(J)
    ranges = local_intertwiner_ranges(J)
    any(isempty, ranges) && return ()

    return collect(Iterators.product(ranges...))
end

function compare_local_vertex(J; max_tuples = nothing)
    tuples = local_intertwiner_tuples(J)

    if max_tuples !== nothing && length(tuples) > max_tuples
        tuples = tuples[1:max_tuples]
    end

    for k_tau in tuples
        value_3j = local_vertex_amplitude(J, k_tau)
        value_6j = local_vertex_amplitude_6j(J, k_tau)

        if abs(value_3j - value_6j) > validation_tolerance
            println("Local 15j mismatch")
            println("J      = ", J)
            println("k_tau  = ", k_tau)
            println("3j     = ", value_3j)
            println("6j     = ", value_6j)
            println("diff   = ", value_6j - value_3j)
            return false
        end
    end

    return true
end

function spin_data_from_local_J(simplex, J)
    spin_keys = oriented_simplex_spin_keys(simplex)

    return Dict(spin_keys[a] => J[a] for a in 1:10)
end

function compare_recoupling_coefficients(simplex, J)
    spin_data = spin_data_from_local_J(simplex, J)
    tetrahedra = oriented_simplex_tetrahedra(simplex)

    for tau in tetrahedra
        global_range = allowed_global_intertwiners(tau, spin_data)
        local_range = allowed_oriented_local_intertwiners(simplex, tau, spin_data)

        for i_tau in global_range, k_tau in local_range
            value_3j = recoupling_R(simplex, tau, spin_data, i_tau, k_tau)
            value_6j = recoupling_R_6j(simplex, tau, spin_data, i_tau, k_tau)

            if abs(value_3j - value_6j) > validation_tolerance
                println("Recoupling coefficient mismatch")
                println("simplex = ", simplex)
                println("tau     = ", tau)
                println("J       = ", J)
                println("(i,k)   = ", (i_tau, k_tau))
                println("3j      = ", value_3j)
                println("6j      = ", value_6j)
                println("diff    = ", value_6j - value_3j)
                return false
            end
        end
    end

    return true
end

function compare_global_vertex(simplex, J; max_tuples = 30)
    spin_data = spin_data_from_local_J(simplex, J)
    tetrahedra = oriented_simplex_tetrahedra(simplex)
    ranges = [allowed_global_intertwiners(tau, spin_data) for tau in tetrahedra]

    if any(isempty, ranges)
        println("Skipping global test with empty global intertwiner range")
        println("simplex = ", simplex)
        println("J       = ", J)
        return true
    end

    tuples = collect(Iterators.product(ranges...))
    if length(tuples) > max_tuples
        tuples = tuples[1:max_tuples]
    end

    for i_tau in tuples
        global_intertwiner_data = Dict(tetrahedra[a] => i_tau[a] for a in 1:5)

        value_3j = simplex_amplitude_recoupled_oriented(
            simplex,
            spin_data,
            global_intertwiner_data,
        )

        value_6j = simplex_amplitude_recoupled_oriented_6j(
            simplex,
            spin_data,
            global_intertwiner_data,
        )

        if abs(value_3j - value_6j) > validation_tolerance
            println("Global vertex mismatch")
            println("simplex = ", simplex)
            println("J       = ", J)
            println("i_tau   = ", i_tau)
            println("3j      = ", value_3j)
            println("6j      = ", value_6j)
            println("diff    = ", value_6j - value_3j)
            return false
        end
    end

    return true
end

function random_local_validation(; trials = 400)
    spin_values = (0 // 1, 1 // 2, 1 // 1, 3 // 2, 2 // 1)
    checked = 0
    attempts = 0

    while checked < trials && attempts < 40 * trials
        attempts += 1
        J = Tuple(rand(spin_values) for _ in 1:10)
        tuples = local_intertwiner_tuples(J)
        isempty(tuples) && continue

        k_tau = rand(tuples)
        value_3j = local_vertex_amplitude(J, k_tau)
        value_6j = local_vertex_amplitude_6j(J, k_tau)

        if abs(value_3j - value_6j) > validation_tolerance
            println("Random local 15j mismatch")
            println("J      = ", J)
            println("k_tau  = ", k_tau)
            println("3j     = ", value_3j)
            println("6j     = ", value_6j)
            println("diff   = ", value_6j - value_3j)
            return false
        end

        checked += 1
    end

    println("random local checks = ", checked)

    return checked == trials
end

function random_global_validation(; trials = 80)
    spin_values = (0 // 1, 1 // 2, 1 // 1, 3 // 2)
    simplices = ((1, 2, 3, 4, 5), (2, 1, 3, 4, 5), (1, 3, 2, 4, 5))
    checked = 0
    attempts = 0

    while checked < trials && attempts < 60 * trials
        attempts += 1
        simplex = rand(simplices)
        J = Tuple(rand(spin_values) for _ in 1:10)
        spin_data = spin_data_from_local_J(simplex, J)
        tetrahedra = oriented_simplex_tetrahedra(simplex)
        ranges = [allowed_global_intertwiners(tau, spin_data) for tau in tetrahedra]

        any(isempty, ranges) && continue

        i_tau = Tuple(rand(range) for range in ranges)
        global_intertwiner_data = Dict(tetrahedra[a] => i_tau[a] for a in 1:5)

        value_3j = simplex_amplitude_recoupled_oriented(
            simplex,
            spin_data,
            global_intertwiner_data,
        )

        value_6j = simplex_amplitude_recoupled_oriented_6j(
            simplex,
            spin_data,
            global_intertwiner_data,
        )

        if abs(value_3j - value_6j) > validation_tolerance
            println("Random global vertex mismatch")
            println("simplex = ", simplex)
            println("J       = ", J)
            println("i_tau   = ", i_tau)
            println("3j      = ", value_3j)
            println("6j      = ", value_6j)
            println("diff    = ", value_6j - value_3j)
            return false
        end

        checked += 1
    end

    println("random global checks = ", checked)

    return checked == trials
end

function main()
    Random.seed!(2026)

    local_test_data = [
        ntuple(_ -> 0 // 1, 10),
        ntuple(_ -> 1 // 2, 10),
        ntuple(_ -> 1 // 1, 10),
        (1 // 1, 1 // 2, 1 // 2, 0 // 1, 1 // 2,
         1 // 2, 0 // 1, 0 // 1, 0 // 1, 0 // 1),
        (1 // 1, 1 // 1, 1 // 2, 1 // 2, 1 // 1,
         1 // 2, 1 // 2, 1 // 2, 1 // 2, 1 // 2),
    ]

    println("Checking fast 6j local vertex against magnetic 3j contraction...")

    local_ok = true
    for J in local_test_data
        ok = compare_local_vertex(J; max_tuples = 250)
        println("local J = ", J, "  ->  ", ok ? "passed" : "failed")
        local_ok &= ok
        local_ok || break
    end

    println()
    println("Checking fast 6j recoupling coefficients...")

    recoupling_ok = true
    for (simplex, J) in [
        ((1, 2, 3, 4, 5), local_test_data[2]),
        ((2, 1, 3, 4, 5), local_test_data[2]),
        ((1, 3, 2, 4, 5), local_test_data[2]),
        ((1, 2, 3, 4, 5), local_test_data[5]),
    ]
        ok = compare_recoupling_coefficients(simplex, J)
        println("R simplex = ", simplex, "  J = ", J, "  ->  ", ok ? "passed" : "failed")
        recoupling_ok &= ok
        recoupling_ok || break
    end

    println()
    println("Checking fast 6j global recoupled oriented vertex...")

    global_ok = true
    global_tests = [
        ((1, 2, 3, 4, 5), local_test_data[2]),
        ((2, 1, 3, 4, 5), local_test_data[2]),
        ((1, 2, 3, 4, 5), local_test_data[5]),
    ]

    for (simplex, J) in global_tests
        ok = compare_global_vertex(simplex, J; max_tuples = 30)
        println("global simplex = ", simplex, "  J = ", J, "  ->  ", ok ? "passed" : "failed")
        global_ok &= ok
        global_ok || break
    end

    println()
    println("Checking random small-spin cases...")

    local_ok &= random_local_validation()
    global_ok &= random_global_validation()

    println()
    println("fast 6j validation summary")
    println("local vertex  : ", local_ok ? "passed" : "failed")
    println("recoupling R  : ", recoupling_ok ? "passed" : "failed")
    println("global vertex : ", global_ok ? "passed" : "failed")

    if !(local_ok && recoupling_ok && global_ok)
        error("fast 6j vertex validation failed")
    end
end

main()
