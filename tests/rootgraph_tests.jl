using Pkg; Pkg.activate("./scripts")
using RootUntangling, RootUntangling.Plots
import RootUntangling: Segment, BranchPoint, PreGraph, RootVertex, RootEdge,
    neighbor, getvertex, roommates, srcs, dsts, cosine_similarity, inner_vertices, outer_vertices

# define root graph
pg = PreGraph(
    [-1, 1, 2, 3, 4],
    [
        [Segment(i, [-1, i]) for i in 1:4];
        Segment(5, [1, 2], 1.0, missing);
        Segment(6, [2, 3], 1.0, missing);
        Segment(7, [2, 4], 1.0, missing);
    ],
    Dict(
        [
            -1 => BranchPoint(-1, :appearance),
            [i => BranchPoint(i, [:a, :b, :c, :d][i], [-1.0, 0.0, 1.0, 0.0][i], [0.0, 0.0, 0.0, 1.0][i]) for i in 1:4]...,
        ]
    )
);
pₛ = 0.1;
rg = get_rootgraph(pg; pₛ, nₕ_min = 1);

# plotting
plot(rg) isa Plots.Plot

# # vertices

# are ids correct (each branchpoint gets 2 root hypotheses at this pₛ)
id.(V₀(rg)) == 1:8
id.(V₊(rg)) == [-1]

# do vertices at the same branchpoint know each other
issetequal(id.(roommates(rg, getvertex(rg, 1))), [1, 2])
issetequal(id.(roommates(rg, getvertex(rg, 3))), [3, 4])

# are neighbors correct
rv = getvertex(rg, 1); # first hypothesis at branchpoint :a
neighbors(rg, rv) isa Vector{<:RootVertex}
issetequal(
    id.(neighbors(rg, rv)),
    [-1, 3, 4]
)

rv = only(V₊(rg));
issetequal(
    id.(neighbors(rg, rv)),
    1:8
)

# inner/outer vertices (only branchpoint :b, holding vertices 3 and 4, has multiple neighboring branchpoints)
issetequal(id.(inner_vertices(rg)), [3, 4])
issetequal(id.(outer_vertices(rg)), [1, 2, 5, 6, 7, 8])

# # edges

# are edges sorted
for edge_set in [E(rg), E₀(rg), E₊(rg)]
    issorted(edge_set, by = src) |> println
end

# are edges unique
for edge_set in [E(rg), E₀(rg), E₊(rg)]
    allunique(edge_set) |> println
end

# are edge sets disjoint
isdisjoint(E₀(rg), E₊(rg))

# do all vertices have a unique set of edges
[allunique(edges(rv)) for rv in V(rg)] |> all

# # segments

# do parallel edges map to the correct segments
length(segments₀(rg)) == 3
all(length.(segments₀(rg)) .== 4) # 2 x 2 root hypotheses per original segment
length(segments₊(rg)) == 4
all(length.(segments₊(rg)) .== 2) # augmented vertex to 2 root hypotheses

seg_ab = only([seg for seg in segments₀(rg) if issetequal(vertices(seg), [1, 2, 3, 4])]);
issetequal(srcs(seg_ab), [1, 2])
issetequal(dsts(seg_ab), [3, 4])
width(seg_ab) == 1.0

# is angle similarity correct
e_ab = first(seg_ab); # edge (1, 3)
seg_bc = only([seg for seg in segments₀(rg) if issetequal(vertices(seg), [3, 4, 5, 6])]);
seg_bd = only([seg for seg in segments₀(rg) if issetequal(vertices(seg), [3, 4, 7, 8])]);
e_bc = first(seg_bc); # edge (3, 5)
e_bd = first(seg_bd); # edge (3, 7)

cosine_similarity(rg, e_ab, e_ab, 1) == 1 # same edge
cosine_similarity(rg, e_ab, e_bc, 3) == -1 # edges with a straight angle
isapprox(cosine_similarity(rg, e_ab, e_bd, 3), 0.0, atol = 1.0e-12) # edges with a 90 degree angle

# # connections
connections = E₂(rg);
allunique(connections)

rv = getvertex(rg, 1);
issetequal(
    [
        vertices.(c) for c in E₂(rv)
    ] .|> sort,
    [
        [(-1, 1), (1, 3)],
        [(-1, 1), (1, 4)],
        [(1, 3), (1, 4)],
    ] .|> sort
)

# # clustering
subgraphs = get_subgraphs(rg; pₛ);
length(subgraphs) == 1
length(V₀(only(subgraphs))) == 8
length(V₊(only(subgraphs))) == 3 # augmentation adds appearance, disappearance and splitting vertices

# # solving (end-to-end on a simple vertical root)
using HiGHS

vertex_data_dict = Dict(
    :a => Dict(:x => 0.0, :y => 2.0, :pred_split => 0.1, :segment_ids => [1]),
    :b => Dict(:x => 0.0, :y => 1.0, :pred_split => 0.1, :segment_ids => [1, 2]),
    :c => Dict(:x => 0.0, :y => 0.0, :pred_split => 0.1, :segment_ids => [2]),
);
edge_data_dict = Dict(
    1 => Dict(:width => 2.0, :pred_primary => 0.9, :xs => [0.0, 0.0], :ys => [2.0, 1.0]),
    2 => Dict(:width => 2.0, :pred_primary => 0.9, :xs => [0.0, 0.0], :ys => [1.0, 0.0]),
);

rg_solve = get_rootgraph(get_pregraph(edge_data_dict, vertex_data_dict; dist_threshold = 0.5); pₛ = 0.5, nₕ_min = 1);
model = solve_rsa(rg_solve; optimizer = HiGHS.Optimizer, num_roots = 1);
roots = get_roots(rg_solve, model);
length(roots) == 1
RootUntangling.is_primary(only(roots))
length(only(roots)) == 3
tortuosity(only(roots)) == 1.0
plot(roots) isa Plots.Plot
plot(rg_solve, get_segment_classification_dict(rg_solve, model)) isa Plots.Plot
hypothesis_plot(rg_solve) isa Plots.Plot
