module tests
using Test, RootUntangling

@testset "PreGraph" include("pregraph_tests.jl")
@testset "SuperGraph" include("supergraph_tests.jl")

end