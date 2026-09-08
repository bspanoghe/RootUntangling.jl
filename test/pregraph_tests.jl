import RootUntangling: segments, getmetavertex, getmetavertices, isspecial, PreGraph, Segment, MetaVertex

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

@test issetequal(vertices(pg), [-3:-1; 1:(length(vertices(pg)) - 3)])
@test length(segments(pg, 1)) == 3
@test length(segments(pg)) == 2 + 2*3 + 1 # 2 real segments, -1 and -2 connect to all real vertices, -3 connects to inner vertex
@test length(getmetavertices(pg)) == 6
@test length([mv for mv in getmetavertices(pg) if isspecial(mv)]) == 3
@test issetequal(neighbors(pg, 1), [-1, -2, 2])
@test issetequal(neighbors(pg, 2), [-1, -2, -3, 1, 3])
@test issetequal(neighbors(pg, 3), [-1, -2, 2])

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
            [i => MetaVertex(i, [:a, :b, :c, :d][i], float(i), float(i)) for i in 1:4]...,
        ]
    )
)

@test issetequal(
    pg.neighbordict,
    Dict(
        [
            -1 => [1, 2, 3, 4],
            1 => [-1, 2],
            2 => [-1, 1, 3, 4],
            3 => [-1, 2],
            4 => [-1, 2],
        ]
    )
)