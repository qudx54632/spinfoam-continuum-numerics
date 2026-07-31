"""BE gauge-volume factor `sum_{l <= j_max} (2l + 1)^2`."""
function one_four_be_factor(j_max)
    return sum(edge_amplitude(l)^2 for l in half_integer_spins(j_max))
end

"""Compute all admissible after-move terms up to one maximum cutoff."""
function one_four_after_terms(boundary_spins, j_max; chi0 = 0)
    return [
        (
            internal_spins = internal_spins,
            needed_cutoff = maximum(internal_spins),
            amplitude = Float64(one_four_after_amplitude(boundary_spins, internal_spins; chi0 = chi0)),
        )
        for internal_spins in one_four_admissible_internal_spins(boundary_spins, j_max)
    ]
end

"""Compare the cutoff sum with `BE factor * before amplitude`."""
function one_four_cutoff_comparison_data(boundary_spins, cutoffs; chi0 = 0)
    max_cutoff = maximum(cutoffs)
    terms = one_four_after_terms(boundary_spins, max_cutoff; chi0 = chi0)
    before = Float64(one_four_before_amplitude(boundary_spins))

    return [
        (
            j_max = Float64(j_max),
            after_sum = sum(
                term.amplitude
                for term in terms
                if term.needed_cutoff <= j_max
            ),
            be_prediction = Float64(one_four_be_factor(j_max)) * before,
        )
        for j_max in cutoffs
    ]
end
