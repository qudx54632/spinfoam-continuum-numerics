const simplex_amplitude_recoupled_oriented_6j_cache = Dict()
const simplex_amplitude_table_recoupled_oriented_6j_cache = Dict()

empty!(simplex_amplitude_recoupled_oriented_6j_cache)
empty!(simplex_amplitude_table_recoupled_oriented_6j_cache)

include(joinpath(@__DIR__, "recoupling_6j.jl"))

"""
Orientation-aware BF 4-simplex vertex in the fixed global tetrahedron basis,
using the fast local 6j expression.

This has the same input/output convention as
`simplex_amplitude_recoupled_oriented`; only the local 15j evaluation is
changed.
"""
function simplex_amplitude_recoupled_oriented_6j(
    simplex,
    spin_data::AbstractDict,
    global_intertwiner_data::AbstractDict,
)
    sigma = oriented_four_simplex(simplex)
    tetrahedra = oriented_simplex_tetrahedra(sigma)

    J = oriented_simplex_spins(sigma, spin_data)
    i_tau = Tuple(global_intertwiner_data[tau] for tau in tetrahedra)

    cache_key = (sigma, J, i_tau)
    if haskey(simplex_amplitude_recoupled_oriented_6j_cache, cache_key)
        return simplex_amplitude_recoupled_oriented_6j_cache[cache_key]
    end

    k_ranges = [
        allowed_oriented_local_intertwiners(sigma, tau, spin_data)
        for tau in tetrahedra
    ]

    if any(isempty, k_ranges)
        simplex_amplitude_recoupled_oriented_6j_cache[cache_key] = 0.0
        return 0.0
    end

    amplitude = 0.0

    for k_tau in Iterators.product(k_ranges...)
        R = recoupling_product_6j(sigma, tetrahedra, spin_data, i_tau, k_tau)
        R == 0.0 && continue

        amplitude += R * local_vertex_amplitude_6j(J, k_tau)
    end

    simplex_amplitude_recoupled_oriented_6j_cache[cache_key] = amplitude

    return amplitude
end

function transform_tensor_axis(tensor, matrix, axis)
    old_dims = size(tensor)
    new_dims = collect(old_dims)
    new_dims[axis] = size(matrix, 1)
    transformed = zeros(Float64, Tuple(new_dims))

    for new_index in CartesianIndices(transformed)
        value = 0.0

        for old_axis_index in 1:size(matrix, 2)
            old_index = ntuple(
                a -> a == axis ? old_axis_index : new_index[a],
                ndims(tensor),
            )
            value += matrix[new_index[axis], old_axis_index] * tensor[old_index...]
        end

        transformed[new_index] = value
    end

    return transformed
end

"""
Compute all recoupled oriented BF amplitudes for a fixed set of global
intertwiner ranges.

Mathematically this evaluates the tensor transform

    A(i_1,...,i_5) = sum_{k_1,...,k_5}
        {15j}_loc(k_1,...,k_5) prod_a R_a(i_a,k_a),

for all requested global intertwiners at once.  It is much faster than calling
`simplex_amplitude_recoupled_oriented_6j` independently for every tuple.
"""
function simplex_amplitude_table_recoupled_oriented_6j(
    simplex,
    spin_data::AbstractDict,
    global_ranges,
)
    sigma = oriented_four_simplex(simplex)
    tetrahedra = oriented_simplex_tetrahedra(sigma)

    J = oriented_simplex_spins(sigma, spin_data)
    global_ranges = Tuple(Tuple(range) for range in global_ranges)

    cache_key = (sigma, J, global_ranges)
    if haskey(simplex_amplitude_table_recoupled_oriented_6j_cache, cache_key)
        return simplex_amplitude_table_recoupled_oriented_6j_cache[cache_key]
    end

    k_ranges = Tuple(
        Tuple(allowed_oriented_local_intertwiners(sigma, tau, spin_data))
        for tau in tetrahedra
    )

    if any(isempty, k_ranges) || any(isempty, global_ranges)
        table = Dict()
        simplex_amplitude_table_recoupled_oriented_6j_cache[cache_key] = table
        return table
    end

    local_tensor = zeros(Float64, Tuple(length(range) for range in k_ranges))

    for k_tuple in Iterators.product(k_ranges...)
        local_index = ntuple(
            a -> findfirst(==(k_tuple[a]), k_ranges[a]),
            length(k_ranges),
        )
        local_tensor[local_index...] = local_vertex_amplitude_6j(J, k_tuple)
    end

    transformed = local_tensor

    for a in 1:5
        matrix = zeros(Float64, length(global_ranges[a]), length(k_ranges[a]))

        for (i_index, i) in enumerate(global_ranges[a]),
            (k_index, k) in enumerate(k_ranges[a])
            matrix[i_index, k_index] = recoupling_R_6j(
                sigma,
                tetrahedra[a],
                spin_data,
                i,
                k,
            )
        end

        transformed = transform_tensor_axis(transformed, matrix, a)
    end

    table = Dict()

    for global_tuple in Iterators.product(global_ranges...)
        global_index = ntuple(
            a -> findfirst(==(global_tuple[a]), global_ranges[a]),
            length(global_ranges),
        )
        table[global_tuple] = transformed[global_index...]
    end

    simplex_amplitude_table_recoupled_oriented_6j_cache[cache_key] = table

    return table
end

"""Fixed-label BF state-sum summand using recoupled oriented 6j vertices."""
function state_sum_amplitude_recoupled_oriented_6j(
    simplices,
    spin_data::AbstractDict,
    global_intertwiner_data::AbstractDict,
    internal_triangles,
)
    amplitude = triangle_amplitude_product(internal_triangles, spin_data)

    for simplex in simplices
        amplitude *= simplex_amplitude_recoupled_oriented_6j(
            simplex,
            spin_data,
            global_intertwiner_data,
        )
    end

    return amplitude
end
