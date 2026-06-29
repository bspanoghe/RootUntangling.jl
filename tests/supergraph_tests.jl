using Pkg; Pkg.activate("./scripts")
using RootUntangling
import RootUntangling: Segment, MetaVertex, PreGraph, Vₕ₀, Vₕ₊, V₀, V₊, V, Vₕ, exclusion_sets, gethypervertex, HyperVertex, SingularVertex, neighbor, getsingularvertex, cosine_similarity, order

# define supergraph
pg = PreGraph(
    [-1, 1, 2, 3, 4],
    [
        [Segment(i, [-1, i]) for i in 1:4];
        Segment(5, [1, 2], 1.0, false);
        Segment(6, [2, 3], 1.0, false);
        Segment(7, [2, 4], 1.0, false);
    ],
    Dict([
        -1 => MetaVertex(-1, :appearance),
        [i => MetaVertex(i, [:a, :b, :c, :d][i], [-1.0, 0.0, 1.0, 0.0][i], [0.0, 0.0, 0.0, 1.0][i]) for i in 1:4]...
    ])
);
n_h = 2;
sg = get_supergraph(pg, n_h);

# plotting
plot(sg) isa Plots.Plot

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

# are exclusion sets correct
hv = Vₕ₀(sg)[1];
issetequal(exclusion_sets(hv), [[1, 2], [1, 3]])

hv = Vₕ₀(sg)[2];
issetequal(exclusion_sets(hv), [[4, 5], [4, 6]])

n_h2 = 3;
sg2 = get_supergraph(pg, n_h2);

hv = Vₕ₀(sg2)[1];
issetequal(exclusion_sets(hv), [[1, 2, 4], [1, 3, 5], [1, 3, 6]])

# gethypervertex 
id(gethypervertex(1, Vₕ₀(sg))) == 1

# order of vertices
hv = Vₕ(sg)[3]
order(hv, vertices(hv)[1]) == 1
order(hv, vertices(hv)[2]) == 2
order(hv, vertices(hv)[3]) == 2

hv = Vₕ(sg2)[3]
order(hv, vertices(hv)[1]) == 1
order(hv, vertices(hv)[2]) == 2
order(hv, vertices(hv)[3]) == 2
order(hv, vertices(hv)[4]) == 3
order(hv, vertices(hv)[5]) == 3
order(hv, vertices(hv)[6]) == 3

# # singular vertices

# are ids correct
id.(V₀(sg)) == 1:length(V₀(sg))
id.(V₊(sg)) == [-1]

# are neighbors correct
sv = V₀(sg)[1];
neighbors(sg, sv) isa Vector{<:SingularVertex}
issetequal(
    id.(neighbors(sg, sv)),
    [-1, 4, 5, 6]
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
vₕ₁, vₕ₂ = Vₕ₀(sg)[[1, 2]];
he = Eₕ(sg)[findfirst(he -> issetequal(vertices(he), (id(vₕ₁), id(vₕ₂))), Eₕ(sg))];
es = E(sg, he);
vs = [getsingularvertex(sg, v) for e in es for v in vertices(e)];
all(hypervertex.(vs) .∈ [[vₕ₁, vₕ₂]])

# do all vertices have a unique set of edges
[allunique(E(v)) for v in V(sg)] |> all

# is angle similarity correct
hes = Eₕ(sg)[[5, 5]] # same edge
cosine_similarity(sg, hes..., 1) == 1

hes = Eₕ(sg)[[5, 6]] # edges with a straight angle
cosine_similarity(sg, hes..., 2) == -1

hes = Eₕ(sg)[[5, 7]] # edges with a 90 degree angle
isapprox(cosine_similarity(sg, hes..., 2), 0.0, atol = 1e-12)

# # connections
connections = E₂(sg);
allunique(connections)

v = V₀(sg)[1];
issetequal(
    [
        vertices.(c) for c in E₂(v)
    ] .|> sort,
    [
        [(-1, 1), (1, 4)],
        [(-1, 1), (1, 5)],
        [(-1, 1), (1, 6)],
        [(1, 4), (1, 5)],
        [(1, 4), (1, 6)],
        [(1, 5), (1, 6)]
    ] .|> sort
)