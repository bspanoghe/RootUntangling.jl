vertex_data_dict = Dict(
    :v1 => Dict(:x => -2, :y => 1, :segment_ids => [1], :pred_split => 0.0),
    :v2 => Dict(:x => -1, :y => 0, :segment_ids => [1, 2, 3], :pred_split => 0.0),
    :v3 => Dict(:x => -2, :y => -1, :segment_ids => [2, 4, 5], :pred_split => 0.0),
    :v4 => Dict(:x => 1, :y => 0, :segment_ids => [3, 6, 7], :pred_split => 0.0),
    :v5 => Dict(:x => -3, :y => -2, :segment_ids => [4, 8, 9], :pred_split => 0.0),
    :v6 => Dict(:x => -2, :y => -2, :segment_ids => [5, 10, 11], :pred_split => 0.0),
    :v7 => Dict(:x => 2, :y => 1, :segment_ids => [7], :pred_split => 0.0),
    :v8 => Dict(:x => 2, :y => -1, :segment_ids => [6, 12, 13], :pred_split => 0.0),
    :v9 => Dict(:x => -4, :y => -3, :segment_ids => [8], :pred_split => 0.0),
    :v10 => Dict(:x => -3, :y => -3, :segment_ids => [9], :pred_split => 0.0),
    :v11 => Dict(:x => -2, :y => -3, :segment_ids => [10], :pred_split => 0.0),
    :v12 => Dict(:x => -1, :y => -3, :segment_ids => [11], :pred_split => 0.0),
    :v13 => Dict(:x => 2, :y => -2, :segment_ids => [12], :pred_split => 0.0),
    :v14 => Dict(:x => 3, :y => -2, :segment_ids => [13], :pred_split => 0.0),
);

edge_data_dict = Dict(
    [i => Dict(:width => 1.0, :pred_primary => 0.5, :xs => Int64[], :ys => Int64[]) for i in 1:13]
);

pg = get_pregraph(edge_data_dict, vertex_data_dict, dist_threshold = 0.5);
rg = get_rootgraph(pg);
model = solve_rsa(rg; optimizer = HiGHS.Optimizer, num_roots = 3, max_roots = 4, ρₒ = 0.3);

# edges correct?
appearance_idxs = src.(E(rg)) .== -1;
@test round(sum(value.(model[:ep₊])[appearance_idxs])) == 3
@test round(sum(value.(model[:ep₋])[appearance_idxs])) == 0

disappearance_idxs = src.(E(rg)) .== -2;
@test round(sum(value.(model[:ep₊])[disappearance_idxs])) == 0
@test round(sum(value.(model[:ep₋])[disappearance_idxs])) == 3

division_idxs = src.(E(rg)) .== -3;
@test round(sum(value.(model[:ep₊])[division_idxs])) == 0
@test round(sum(value.(model[:ep₋])[division_idxs])) == 0

# connections correct?
appearance_idxs = in.([-1], vertices.(first.(E₂(rg))));
@test round(sum(value.(model[:cp₊])[appearance_idxs])) == 3
@test round(sum(value.(model[:cp₋])[appearance_idxs])) == 0

disappearance_idxs = in.([-2], vertices.(first.(E₂(rg))));
@test round(sum(value.(model[:cp₊])[disappearance_idxs])) == 0
@test round(sum(value.(model[:cp₋])[disappearance_idxs])) == 3

division_idxs = in.([-3], vertices.(first.(E₂(rg))));
@test round(sum(value.(model[:cp₊])[division_idxs])) == 0
@test round(sum(value.(model[:cp₋])[division_idxs])) == 0

# roots
rootsystems = get_rootsystems(rg, model);
@test length(rootsystems) == 3

rootsystems_new = greedy_switch(rg, rootsystems);
@test length(rootsystems) == length(rootsystems_new)
@test length(reduce(vcat, laterals.(rootsystems))) == length(reduce(vcat, laterals.(rootsystems_new)))
@test sum(curve_length.([rg], reduce(vcat, laterals.(rootsystems)))) ==
    sum(curve_length.([rg], reduce(vcat, laterals.(rootsystems_new))))
