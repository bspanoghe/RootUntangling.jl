vertex_data_dict = Dict(
    :a1 => Dict(:x => -1, :y => -1, :segment_ids => [1, 3, 4], :pred_split => 0.0),
    :a2 => Dict(:x => -2, :y => -2, :segment_ids => [4, 7, 8], :pred_split => 0.0),
    :a3 => Dict(:x => -3, :y => -2, :segment_ids => [7], :pred_split => 0.0),
    :a4 => Dict(:x => -2, :y => -3, :segment_ids => [8], :pred_split => 0.0),
    :b1 => Dict(:x => 0, :y => 1, :segment_ids => [1, 2, 5], :pred_split => 0.0),
    :b2 => Dict(:x => 0, :y => 2, :segment_ids => [5, 9, 10], :pred_split => 0.0),
    :b3 => Dict(:x => -0.5, :y => 3, :segment_ids => [9], :pred_split => 0.0),
    :b4 => Dict(:x => 0.5, :y => 3, :segment_ids => [10], :pred_split => 0.0),
    :c1 => Dict(:x => 1, :y => -1, :segment_ids => [2, 3, 6], :pred_split => 0.0),
    :c2 => Dict(:x => 2, :y => -2, :segment_ids => [6, 11, 12], :pred_split => 0.0),
    :c3 => Dict(:x => 3, :y => -2, :segment_ids => [11], :pred_split => 0.0),
    :c4 => Dict(:x => 2, :y => -3, :segment_ids => [12], :pred_split => 0.0),
);

edge_data_dict = Dict(
    [i => Dict(:width => 1.0, :pred_primary => 0.5, :xs => Int64[], :ys => Int64[]) for i in 1:12]
);

pg = get_pregraph(edge_data_dict, vertex_data_dict, dist_threshold = 0.5);
rg = get_rootgraph(pg);
graphplot(rg)

model = solve_rsa(rg; optimizer = HiGHS.Optimizer, time_limit = 60, num_roots = 6,
    w_gp = 0.0, w_gl = 0.0, max_roots = 2, ρₘ_max = 0.99
)
graphplot(rg, model, augmented_alpha = 0.3, size = (400, 300), vertex_kwargs = Dict(:markersize => 2))

roots = get_rootsystems(rg, model);
rootplot(rg, roots)

rootplot(rg, roots[1])
rootplot(rg, roots[2])
rootplot(rg, roots[3])
rootplot(rg, roots[4])
rootplot(rg, roots[5])
rootplot(rg, roots[6])