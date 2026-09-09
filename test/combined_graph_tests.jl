vertex_data_dict = Dict(
    :a => Dict(:x => -1, :y => 0, :segment_ids => [1], :pred_split => 0.0),
    :b => Dict(:x => 0, :y => 1, :segment_ids => [1, 2], :pred_split => 0.0),
    :c => Dict(:x => 1, :y => 0, :segment_ids => [2], :pred_split => 0.0),
    :c2 => Dict(:x => 1, :y => 0, :segment_ids => [2], :pred_split => 0.0),
    :isolated => Dict(:x => 1, :y => 1, :segment_ids => [3], :pred_split => 0.0)
)

edge_data_dict = Dict(
    1 => Dict(:width => 1.5, :pred_primary => 0.0, :xs => [1, 2], :ys => [1, 2]),
    2 => Dict(:width => 2.0, :pred_primary => 0.0, :xs => [1, 2], :ys => [1, 2]),
    3 => Dict(:width => 0.5, :pred_primary => 0.0, :xs => [1, 2], :ys => [1, 2]),
)

pg = get_pregraph(edge_data_dict, vertex_data_dict, dist_threshold = 0.5);
rg = get_supergraph(pg);

# # vertices
@test length(V₀(rg)) == 3
@test length(V₊(rg)) == 3
@test issetequal(
    V(rg),
    [V₀(rg); V₊(rg)]
)

@test id.(V₀(rg)) == 1:3
@test id.(V₊(rg)) == -1:-1:-3

# # edges
# ## direction
rv = V₊(rg)[1]
@test all([direction(rv, e) == 1 for e in edges(rv)])

rv = V₊(rg)[2]
@test all([direction(rv, e) == 1 for e in edges(rv)])

rvv = V₊(rg)[3]
@test all([direction(rv, e) == 1 for e in edges(rv)])