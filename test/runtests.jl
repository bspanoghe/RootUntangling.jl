module tests
using Test, RootUntangling, JuMP, HiGHS

@testset "Pregraph" include("pregraph_tests.jl")
@testset "Supergraph" include("rootgraph_tests.jl")
@testset "Combined graphs" include("combined_graph_tests.jl")
@testset "Solving" include("solving_tests.jl")

end