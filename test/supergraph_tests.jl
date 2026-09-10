import RootUntangling: Segment, MetaVertex, PreGraph, V₀, V₊, V,
    RootVertex, RootEdge, neighbor, getrootvertex,
    cosine_similarity, direction
import RootUntangling.Makie: Figure

# define rootgraph
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
rg = get_rootgraph(pg);

# plotting
@test graphplot(rg) isa Figure

# # vertices
# are ids correct
@test issetequal(id.(V₀(rg)), 1:4)
@test id(only(V₊(rg))) == -1

# are neighbors correct
rv = V₀(rg)[1];
@test neighbors(rg, rv) isa Vector{<:RootVertex}
@test issetequal(
    id.(neighbors(rg, rv)),
    [-1, 2]
)

rv = only(V₊(rg));
@test issetequal(
    id.(neighbors(rg, rv)),
    [1, 2, 3, 4]
)

# # edges
# correct direction
rv = V₀(rg)[1]

e = edges(rv)[findfirst(is_augmented, edges(rv))]
@test direction(rv, e) == -1

e = edges(rv)[findfirst(e -> e isa RootEdge && !is_augmented(e), edges(rv))]
@test direction(rv, e) == 1

es = filter(e -> e isa RootEdge, edges(rv))
@test issetequal(
    [direction(rv, e) for e in es],
    [1, -1]
)

rv = V₊(rg)[1]
@test all([direction(rv, e) == 1 for e in edges(rv)])

# are edges sorted
for edge_set in [E(rg), E₀(rg), E₊(rg)]
    @test issorted(edge_set, by = e -> src(e))
end

# are edges unique
for edge_set in [E(rg), E₀(rg), E₊(rg)]
    @test allunique(edge_set)
end

# are edge sets disjoint
@test isdisjoint(E₀(rg), E₊(rg))

# do all vertices have a unique set of edges
@test [allunique(E(v)) for v in V(rg)] |> all

# is angle similarity correct
res = E(rg)[[5, 5]] # same edge
@test cosine_similarity(rg, res..., 1) == 1

res = E(rg)[[5, 6]] # edges with a straight angle
@test cosine_similarity(rg, res..., 2) == -1

res = E(rg)[[5, 7]] # edges with a 90 degree angle
@test isapprox(cosine_similarity(rg, res..., 2), 0.0, atol = 1.0e-12)

# # connections
connections = E₂(rg);
@test allunique(connections)

issetequal(
    unique(reduce(vcat, E₂.([rg], E₀(rg)))),
    E₂(rg)
)

rv = V₀(rg)[1];
@test issetequal(
    [
        vertices.(c) for c in E₂(rv)
    ] .|> sort,
    [
        [(-1, 1), (1, 2)],
    ] .|> sort
)