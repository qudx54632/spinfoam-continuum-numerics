include(joinpath(@__DIR__, "..", "core", "load_bf4d.jl"))

const tolerance = 1.0e-10

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

function compare_local_vertex(J)
    for k_tau in local_intertwiner_tuples(J)
        old_value = local_vertex_amplitude(J, k_tau)
        fast_value = local_vertex_amplitude_6j(J, k_tau)

        if abs(old_value - fast_value) > tolerance
            println("Local 15j mismatch")
            println("J      = ", J)
            println("k_tau  = ", k_tau)
            println("old 3j = ", old_value)
            println("fast6j = ", fast_value)
            println("diff   = ", fast_value - old_value)
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
            old_value = recoupling_R(simplex, tau, spin_data, i_tau, k_tau)
            fast_value = recoupling_R_6j(simplex, tau, spin_data, i_tau, k_tau)

            if abs(old_value - fast_value) > tolerance
                println("Recoupling coefficient mismatch")
                println("simplex = ", simplex)
                println("tau     = ", tau)
                println("J       = ", J)
                println("(i,k)   = ", (i_tau, k_tau))
                println("old 3j  = ", old_value)
                println("fast6j  = ", fast_value)
                println("diff    = ", fast_value - old_value)
                return false
            end
        end
    end

    return true
end

function compare_global_vertex(simplex, J)
    spin_data = spin_data_from_local_J(simplex, J)
    tetrahedra = oriented_simplex_tetrahedra(simplex)
    ranges = [allowed_global_intertwiners(tau, spin_data) for tau in tetrahedra]

    any(isempty, ranges) && return true

    for i_tau in Iterators.product(ranges...)
        global_intertwiner_data = Dict(tetrahedra[a] => i_tau[a] for a in 1:5)

        old_value = simplex_amplitude_recoupled_oriented(
            simplex,
            spin_data,
            global_intertwiner_data,
        )

        fast_value = simplex_amplitude_recoupled_oriented_6j(
            simplex,
            spin_data,
            global_intertwiner_data,
        )

        if abs(old_value - fast_value) > tolerance
            println("Global vertex mismatch")
            println("simplex = ", simplex)
            println("J       = ", J)
            println("i_tau   = ", i_tau)
            println("old 3j  = ", old_value)
            println("fast6j  = ", fast_value)
            println("diff    = ", fast_value - old_value)
            return false
        end
    end

    return true
end

function main()
    test_data = [
        ntuple(_ -> 0 // 1, 10),
        ntuple(_ -> 1 // 2, 10),
        ntuple(_ -> 1 // 1, 10),
        (1 // 1, 1 // 2, 1 // 2, 0 // 1, 1 // 2,
         1 // 2, 0 // 1, 0 // 1, 0 // 1, 0 // 1),
        (1 // 1, 1 // 1, 1 // 2, 1 // 2, 1 // 1,
         1 // 2, 1 // 2, 1 // 2, 1 // 2, 1 // 2),
    ]

    println("Checking BF_4d fast 6j local vertex against old 3j vertex...")

    local_ok = true
    for J in test_data
        ok = compare_local_vertex(J)
        println("local J = ", J, "  ->  ", ok ? "passed" : "failed")
        local_ok &= ok
        local_ok || break
    end

    println()
    println("Checking BF_4d fast 6j recoupling coefficients...")

    recoupling_ok = true
    for (simplex, J) in [
        ((1, 2, 3, 4, 5), test_data[2]),
        ((2, 1, 3, 4, 5), test_data[2]),
        ((1, 3, 2, 4, 5), test_data[2]),
        ((1, 2, 3, 4, 5), test_data[5]),
    ]
        ok = compare_recoupling_coefficients(simplex, J)
        println("R simplex = ", simplex, "  J = ", J, "  ->  ", ok ? "passed" : "failed")
        recoupling_ok &= ok
        recoupling_ok || break
    end

    println()
    println("Checking BF_4d fast 6j global oriented recoupled vertex...")

    global_ok = true
    global_tests = [
        ((1, 2, 3, 4, 5), test_data[2]),
        ((2, 1, 3, 4, 5), test_data[2]),
        ((1, 2, 3, 4, 5), test_data[5]),
    ]

    for (simplex, J) in global_tests
        ok = compare_global_vertex(simplex, J)
        println("global simplex = ", simplex, "  J = ", J, "  ->  ", ok ? "passed" : "failed")
        global_ok &= ok
        global_ok || break
    end

    println()
    println("fast 6j validation summary")
    println("local vertex  : ", local_ok ? "passed" : "failed")
    println("recoupling R  : ", recoupling_ok ? "passed" : "failed")
    println("global vertex : ", global_ok ? "passed" : "failed")

    if !(local_ok && recoupling_ok && global_ok)
        error("BF_4d fast 6j validation failed")
    end
end

main()
