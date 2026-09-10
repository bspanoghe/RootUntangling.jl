module tests
using Test, RootUntangling

@testset "Pregraph" include("pregraph_tests.jl")
@testset "Supergraph" include("rootgraph_tests.jl")
@testset "Combined graphs" include("combined_graph_tests.jl")

end