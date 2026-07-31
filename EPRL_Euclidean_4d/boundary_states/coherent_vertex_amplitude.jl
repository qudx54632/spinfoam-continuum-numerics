"""Single vertices coupled to Livine-Speziale coherent intertwiners."""

"""Precompute the coherent boundary coefficient `c_i(tau)` for one vertex."""
function coherent_intertwiner_coefficient_table(
    tetrahedra,
    intertwiner_ranges,
    spin_data::AbstractDict,
    normal_data::AbstractDict,
)
    coefficients = Dict()

    for (tetrahedron, i_values) in zip(tetrahedra, intertwiner_ranges)
        for i in i_values
            coefficients[(tetrahedron, i)] = coherent_intertwiner_coefficient(
                tetrahedron,
                spin_data,
                normal_data,
                i,
            )
        end
    end

    return coefficients
end

"""Read the product `prod_tau c_{i_tau}(tau)` from a precomputed table."""
function coherent_intertwiner_coefficient_product(
    coefficient_table::AbstractDict,
    tetrahedra,
    intertwiners,
)
    product = 1.0 + 0.0im

    for (tetrahedron, i) in zip(tetrahedra, intertwiners)
        product *= coefficient_table[(tetrahedron, i)]
        abs(product) == 0.0 && return product
    end

    return product
end

function coherent_su2_bf_single_vertex_amplitude_with_vertex(
    vertex_amplitude,
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

    total = 0.0 + 0.0im

    for intertwiners in Iterators.product(intertwiner_ranges...)
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

        total += vertex_amplitude(
            sigma,
            spin_data,
            intertwiner_data,
        ) * coefficient_product
    end

    return total
end

"""
Fixed-spin coherent SU(2) BF vertex

    W_v^{SU(2)}(j_f,n_{tau f})
      = sum_{i_tau} A_v^{BF}(j_f,i_tau)
        prod_tau c_{i_tau}(j_{f subset tau}, n_{tau f}).

This is the clean SU(2) 15j object used before adding the EPRL fusion map.
The normal data should preferably be keyed by `(tetrahedron, triangle)`.
"""
function coherent_su2_bf_single_vertex_amplitude(
    simplex,
    spin_data::AbstractDict,
    normal_data::AbstractDict,
)
    return coherent_su2_bf_single_vertex_amplitude_with_vertex(
        simplex_amplitude_recoupled_oriented,
        simplex,
        spin_data,
        normal_data,
    )
end

"""Same coherent SU(2) BF vertex, using the fast 6j local 15j backend."""
function coherent_su2_bf_single_vertex_amplitude_6j(
    simplex,
    spin_data::AbstractDict,
    normal_data::AbstractDict,
)
    return coherent_su2_bf_single_vertex_amplitude_with_vertex(
        simplex_amplitude_recoupled_oriented_6j,
        simplex,
        spin_data,
        normal_data,
    )
end

"""
Fixed-spin coherent EPRL vertex

    W_v(j_f,n_{tau f})
      = sum_{i_tau} A_v(j_f,i_tau)
        prod_tau c_{i_tau}(j_{f subset tau}, n_{tau f}).

The normal data should preferably be keyed by `(tetrahedron, triangle)`.
"""
function coherent_single_vertex_amplitude(
    simplex,
    spin_data::AbstractDict,
    normal_data::AbstractDict,
    gamma,
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

    total = 0.0 + 0.0im

    for intertwiners in Iterators.product(intertwiner_ranges...)
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

        total += eprl_vertex_amplitude(
            sigma,
            spin_data,
            intertwiner_data,
            gamma,
        ) * coefficient_product
    end

    return total
end
