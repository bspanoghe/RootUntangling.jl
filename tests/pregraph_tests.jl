import RootUntangling: segments, getmetavertex, getmetavertices, isspecial, PreGraph, Segment, MetaVertex

vertex_data_dict = Dict(
    :a => Dict(:x => -1, :y => 0, :segment_ids => [1]),
    :b => Dict(:x => 0, :y => 1, :segment_ids => [1, 2]),
    :c => Dict(:x => 1, :y => 0, :segment_ids => [2]),
    :c2 => Dict(:x => 1, :y => 0, :segment_ids => [2]),
    :isolated => Dict(:x => 1, :y => 1, :segment_ids => [3])
)

edge_data_dict = Dict(
    1 => Dict(:fake => false, :width => 1.5),
    2 => Dict(:fake => false, :width => 2.0),
    3 => Dict(:fake => false, :width => 0.5),
)

pg = get_pregraph(edge_data_dict, vertex_data_dict, dist_threshold = 0.5);

issetequal(vertices(pg), [-3:-1; 1:(length(vertices(pg))-3)])
length(segments(pg, 1)) == 4
length(segments(pg)) == 3*3 + 2
length(getmetavertices(pg)) == 6
length([mv for mv in getmetavertices(pg) if isspecial(mv)]) == 3
issetequal(neighbors(pg, 1), [-1, -2, -3, 2])
issetequal(neighbors(pg, 2), [-1, -2, -3, 1, 3])
issetequal(neighbors(pg, 3), [-1, -2, -3, 2])

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
        [i => MetaVertex(i, [:a, :b, :c, :d][i], float(i), float(i)) for i in 1:4]...
    ])
)

issetequal(
    pg.neighbordict,
    Dict([
        -1 => [1, 2, 3, 4],
        1 => [-1, 2],
        2 => [-1, 1, 3, 4],
        3 => [-1, 2],
        4 => [-1, 2]
    ])
)