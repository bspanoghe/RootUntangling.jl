import RootUntangling: Segment, MetaVertex, PreGraph, Vₕ₀, Vₕ₊, V₀, V₊, V, Vₕ,
    gethypervertex, HyperVertex, SingularVertex, ExternalEdge, InternalEdge, neighbor, getsingularvertex,
    cosine_similarity, polarity, direction
import RootUntangling.Makie: Figure

# define supergraph
pg = PreGraph(
    [-1, 1, 2, 3, 4],
    [
        [Segment(i, [-1, i]) for i in 1:4];
        Segment(5, [1, 2]);
        Segment(6, [2, 3]);
        Segment(7, [2, 4]);
    ],
    Dict(
        [
            -1 => MetaVertex(-1, :appearance),
            [i => MetaVertex(i, [:a, :b, :c, :d][i], [-1.0, 0.0, 1.0, 0.0][i], [0.0, 0.0, 0.0, 1.0][i]) for i in 1:4]...,
        ]
    )
);
sg = get_supergraph(pg);

# plotting
@test graphplot(sg) isa Figure

# # hypervertices
# are ids correct
@test issetequal(id.(Vₕ₀(sg)), 1:4)
@test id(only(Vₕ₊(sg))) == -1

# are neighbors correct
hv = Vₕ₀(sg)[1];
@test neighbors(sg, hv) isa Vector{<:HyperVertex}
@test issetequal(
    id.(neighbors(sg, hv)),
    [-1, 2]
)

hv = only(Vₕ₊(sg));
@test issetequal(
    id.(neighbors(sg, hv)),
    [1, 2, 3, 4]
)

# gethypervertex
@test id(gethypervertex(Vₕ₀(sg), 1)) == 1

# # singular vertices

# are ids correct
@test id.(V₀(sg)) == 1:length(V₀(sg))
@test id.(V₊(sg)) == [-1]

# are neighbors correct
sv = V₀(sg)[1];
@test neighbors(sg, sv) isa Vector{<:SingularVertex}
@test issetequal(
    id.(neighbors(sg, sv)),
    [-1, 2, 3]
)
@test allunique(id.(neighbors(sg, sv)))

sv = V(sg, Vₕ₀(sg)[2])[2];
@test neighbors(sg, sv) isa Vector{<:SingularVertex}
@test issetequal(
    id.(neighbors(sg, sv)),
    [-1, 2, 3, 6, 8]
)
@test allunique(id.(neighbors(sg, sv)))

sv = V₊(sg)[1];
@test issetequal(
    neighbors(sg, sv),
    V₀(sg)
)
@test allunique(id.(neighbors(sg, sv)))

# # edges

# are the correct edges in place
sv = V₀(sg)[1]
@test issetequal(
    vertices.(edges(sv)),
    [(-1, 1), (1, 2), (2, 1), (1, 3)]
)
@test allunique(vertices.(edges(sv)))

# correct polarity
sv = V₀(sg)[1]
@test sum(polarity.(edges(sv)) .== -1) == 1 # only intra edge has negative polarity

sv = V₀(sg)[2]
@test sum(polarity.(edges(sv)) .== 1) == 2 # only augmented edge and intra edge have positive polarity

# correct direction
sv = V₀(sg)[1]

e = edges(sv)[findfirst(is_augmented, edges(sv))]
@test direction(sv, e) == -1

e = edges(sv)[findfirst(e -> e isa ExternalEdge && !is_augmented(e), edges(sv))]
@test direction(sv, e) == 1

es = filter(e -> e isa InternalEdge, edges(sv))
@test issetequal(
    [direction(sv, e) for e in es],
    [1, -1]
)



# are edges sorted
for edge_set in [E(sg), E₀(sg), E₊(sg), Eₕ(sg), Eₕ₀(sg), Eₕ₊(sg)]
    @test issorted(edge_set, by = e -> src(e))
end

# are edges unique
for edge_set in [E(sg), E₀(sg), E₊(sg), Eₕ(sg), Eₕ₀(sg), Eₕ₊(sg)]
    @test allunique(edge_set)
end

# are edge sets disjoint
@test isdisjoint(E₀(sg), E₊(sg))
@test isdisjoint(Eₕ₀(sg), Eₕ₊(sg))

# do singular edges map to the correct hyperedges
hv1, hv2 = Vₕ₀(sg)[[1, 2]];
he = Eₕ(sg)[findfirst(he -> issetequal(vertices(he), (id(hv1), id(hv2))), Eₕ(sg))]; # he between hvs
es = E(sg, he); # singular edges of he
svs = [getsingularvertex(sg, v) for e in es for v in vertices(e)]; # singular vertices of ses
@test all(hypervertex.(svs) .∈ [[hv1, hv2]])

# do all vertices have a unique set of edges
@test [allunique(E(v)) for v in V(sg)] |> all

# is angle similarity correct
hes = Eₕ(sg)[[5, 5]] # same edge
@test cosine_similarity(sg, hes..., 1) == 1

hes = Eₕ(sg)[[5, 6]] # edges with a straight angle
@test cosine_similarity(sg, hes..., 2) == -1

hes = Eₕ(sg)[[5, 7]] # edges with a 90 degree angle
@test isapprox(cosine_similarity(sg, hes..., 2), 0.0, atol = 1.0e-12)

# # connections
connections = E₂(sg);
@test allunique(connections)

sv = V₀(sg)[1];
@test issetequal(
    [
        vertices.(c) for c in E₂(sv)
    ] .|> sort,
    [
        [(-1, 1), (1, 2)],
        [(-1, 1), (2, 1)],
        [(-1, 1), (1, 3)],
        [(1, 2), (1, 3)],
        [(2, 1), (1, 3)],
    ] .|> sort
)