function get_supergraph(pg::PreGraph; pₛ = 0.2)
    n_hs = [get_num_hypotheses(pg, mv; pₛ) for mv in getmetavertices(pg) if !isspecial(mv)]

    Vₕ₀ = [
        get_hypervertex(pg, mv, n_hs, i)
        for (i, mv) in enumerate(getmetavertices(pg))
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

# get amount of hypotheses corresponding to a metavertex
function get_num_hypotheses(pg::PreGraph, mv::MetaVertex; pₛ)
    n_h = [get_num_hypotheses(pg, s; pₛ) for s in segments(pg, mv) if !isspecial(s)] |> maximum

    return n_h
end

function get_num_hypotheses(pg::PreGraph, s::Segment; pₛ)
    all_widths = [width(s) for s in segments(pg) if !isspecial(s)]
    single_width = quantile(all_widths, pₛ)
    n_h = div(width(s), single_width, RoundUp) |> Int

    return n_h
end

# instantiate a hypervertex based on a metavertex
function get_hypervertex(pg::PreGraph{T, U, V}, mv::MetaVertex{T, U}, n_hs::Vector{<:Integer}, i::Integer) where {T, U, V}
    prev_id = sum(n_hs[1:i-1]) # amount of vertices that have been defined in previous hypervertices
    vertices = [i for i in prev_id+1:prev_id+n_hs[i]]

    return HyperVertex(id(mv), HyperEdge.(segments(pg, mv)), x(mv), y(mv), pred_split(mv), vertices)
end

HyperEdge(s::Segment) = HyperEdge(vertices(s)..., width(s), pred_primary(s))
SingularEdge(s::Segment) = SingularEdge(vertices(s)..., HyperEdge(s))

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
        pₛ::Real, flip_y::Bool
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
    sg = get_supergraph(pg; pₛ)

    return sg
end