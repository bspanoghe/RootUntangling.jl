function get_supergraph(pg::PreGraph, n_h::Integer)
    Vₕ₀ = [
        HyperVertex(id(mv), HyperEdge.(segments(pg, mv)), x(mv), y(mv), pred_split(mv), n_h)
        for mv in getmetavertices(pg)
        if !isspecial(mv)
    ]
    Vₕ₊ = [
        HyperVertex(id(mv), HyperEdge.(segments(pg, mv)), x(mv), y(mv), NaN, [id(mv)])
        for mv in getmetavertices(pg)
        if isspecial(mv)
    ]
    V₀ = (
    [
        getsingularvertices(hv, [Vₕ₀; Vₕ₊])
        for hv in Vₕ₀
    ] |> x -> reduce(vcat, x)
    )
    V₊ = [
        getsingularvertices(hv, [Vₕ₀; Vₕ₊])
        for hv in Vₕ₊
    ] |> x -> reduce(vcat, x)

    return SuperGraph(Vₕ₀, Vₕ₊, V₀, V₊)
end

HyperEdge(s::Segment) = HyperEdge(vertices(s)..., width(s), pred_primary(s))
SingularEdge(s::Segment) = SingularEdge(vertices(s)..., HyperEdge(s))

# instantiate a hypervertex based on the number of hypotheses `n_h`
function HyperVertex(id::T, edges::Vector{HyperEdge{T, U}}, x::U, y::U, pred_split::Union{U, Missing}, n_h::Integer) where {T, U}
    n_v = sum(1:n_h) # number of vertices per hypervertex
    prev_id = (id-1) * n_v # amount of vertices that have been defined in previous hypervertices
    vertices = [i for i in prev_id+1:prev_id+n_v]

    return HyperVertex(id, edges, x, y, pred_split, vertices)
end

# get singular vertices corresponding with a given hypervertex (containing the correct edges)
function getsingularvertices(hv::HyperVertex{T, U}, Vₕ::Vector{HyperVertex{T, U}}) where {T, U}
    [
        SingularVertex(
            v,
            [
                SingularEdge{T, U}[
                    SingularEdge(v, v_nb, he) 
                    for v_nb in vertices(neighbor(hv, he, Vₕ))
                ]
                for he in edges(hv)
            ] |> x -> reduce(vcat, x, init = SingularEdge{T, U}[]),
            hv
        )
        for v in vertices(hv)
    ]
end

# method for doing all steps from file reading
function get_supergraph(filename_segments::String, filename_vertices::String; dist_threshold::Real,
        num_hypotheses::Real, flip_y::Bool
    )
    
    edge_data = read_data(filename_segments, :Segment_ID)
    edge_data_dict = get_edge_info(edge_data, dist_colname = :Mean_Distance, fake_colname = :Fake_Lateral,
        primary_score_colname = :Heatmap_Mean
    )

    vertex_data = read_data(filename_vertices, :Node)
    vertex_data_dict = get_vertex_info(vertex_data; flip_y, segment_ids_colname = :Segment_IDs,
        x_colname = :Coord_x, y_colname = :Coord_y, lateral_score_colname = :Lateral_Score
    )

    pg = get_pregraph(edge_data_dict, vertex_data_dict; dist_threshold)
    sg = get_supergraph(pg, num_hypotheses)

    return sg
end