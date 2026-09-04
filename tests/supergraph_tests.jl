using Pkg; Pkg.activate("./scripts")
using RootUntangling
import RootUntangling: Segment, MetaVertex, PreGraph, Vₕ₀, Vₕ₊, V₀, V₊, V, Vₕ, gethypervertex, HyperVertex, SingularVertex, neighbor, getsingularvertex, cosine_similarity

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
graphplot(sg) isa Figure

# # hypervertices
# are ids correct
issetequal(id.(Vₕ₀(sg)), 1:4)
id(only(Vₕ₊(sg))) == -1

# are neighbors correct
hv = Vₕ₀(sg)[1];
neighbors(sg, hv) isa Vector{<:HyperVertex}
issetequal(
    id.(neighbors(sg, hv)),
    [-1, 2]
)

hv = only(Vₕ₊(sg));
issetequal(
    id.(neighbors(sg, hv)),
    [1, 2, 3, 4]
)

# gethypervertex
id(gethypervertex(Vₕ₀(sg), 1)) == 1

# # singular vertices

# are ids correct
id.(V₀(sg)) == 1:length(V₀(sg))
id.(V₊(sg)) == [-1]

# are neighbors correct
sv = V₀(sg)[1];
neighbors(sg, sv) isa Vector{<:SingularVertex}
issetequal(
    id.(neighbors(sg, sv)),
    [-1, 2, 3]
)

sv = V(sg, Vₕ₀(sg)[2])[2];
neighbors(sg, sv) isa Vector{<:SingularVertex}
issetequal(
    id.(neighbors(sg, sv)),
    [-1, 2, 3, 6, 8]
)

sv = V₊(sg)[1];
issetequal(
    neighbors(sg, sv),
    V₀(sg)
)

# # edges

# are edges sorted
for edge_set in [E(sg), E₀(sg), E₊(sg), Eₕ(sg), Eₕ₀(sg), Eₕ₊(sg)]
    issorted(edge_set, by = e -> src(e)) |> println
end

# are edges unique
for edge_set in [E(sg), E₀(sg), E₊(sg), Eₕ(sg), Eₕ₀(sg), Eₕ₊(sg)]
    allunique(edge_set) |> println
end

# are edge sets disjoint
isdisjoint(E₀(sg), E₊(sg))
isdisjoint(Eₕ₀(sg), Eₕ₊(sg))

# do singular edges map to the correct hyperedges
hv1, hv2 = Vₕ₀(sg)[[1, 2]];
he = Eₕ(sg)[findfirst(he -> issetequal(vertices(he), (id(hv1), id(hv2))), Eₕ(sg))]; # he between hvs
es = E(sg, he); # singular edges of he
svs = [getsingularvertex(sg, v) for e in es for v in vertices(e)]; # singular vertices of ses
all(hypervertex.(svs) .∈ [[hv1, hv2]])

# do all vertices have a unique set of edges
[allunique(E(v)) for v in V(sg)] |> all

# is angle similarity correct
hes = Eₕ(sg)[[5, 5]] # same edge
cosine_similarity(sg, hes..., 1) == 1

hes = Eₕ(sg)[[5, 6]] # edges with a straight angle
cosine_similarity(sg, hes..., 2) == -1

hes = Eₕ(sg)[[5, 7]] # edges with a 90 degree angle
isapprox(cosine_similarity(sg, hes..., 2), 0.0, atol = 1.0e-12)

# # connections
connections = E₂(sg);
allunique(connections)

sv = V₀(sg)[1];
issetequal(
    [
        vertices.(c) for c in E₂(sv)
    ] .|> sort,
    [
        [(-1, 1), (1, 2)],
        [(-1, 1), (1, 3)],
        [(1, 2), (1, 3)],
    ] .|> sort
)
