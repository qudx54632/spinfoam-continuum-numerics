include(joinpath(@__DIR__, "core", "spin_utils.jl"))
include(joinpath(@__DIR__, "core", "state_sum_utils.jl"))
include(joinpath(@__DIR__, "core", "orientation.jl"))

include(joinpath(@__DIR__, "su2_bf_vertex", "intertwiner_basis.jl"))
include(joinpath(@__DIR__, "su2_bf_vertex", "local_vertex_amplitude.jl"))
include(joinpath(@__DIR__, "su2_bf_vertex", "recoupling.jl"))
include(joinpath(@__DIR__, "su2_bf_vertex", "global_amplitude.jl"))

include(joinpath(@__DIR__, "su2_bf_vertex_fast_6j", "local_vertex_amplitude_6j.jl"))
include(joinpath(@__DIR__, "su2_bf_vertex_fast_6j", "global_amplitude_6j.jl"))

include(joinpath(@__DIR__, "eprl_vertex", "eprl_simplicity.jl"))
include(joinpath(@__DIR__, "eprl_vertex", "fusion_coefficient.jl"))
include(joinpath(@__DIR__, "eprl_vertex", "eprl_vertex_amplitude.jl"))

include(joinpath(@__DIR__, "boundary_states", "coherent_intertwiner.jl"))
include(joinpath(@__DIR__, "boundary_states", "geometric_normals.jl"))
include(joinpath(@__DIR__, "boundary_states", "coherent_vertex_amplitude.jl"))

include(joinpath(@__DIR__, "complexes", "complex_amplitude.jl"))
include(joinpath(@__DIR__, "complexes", "one_five_move.jl"))
