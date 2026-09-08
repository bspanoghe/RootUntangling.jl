"""
    get_supergraph(filename_segments::String, filename_vertices::String;
        dist_threshold::Real, reverse_y::Bool, pₛ::Real = 0.2, nₕ_min::Integer = 1, [colnames...]
    )

Create a supergraph from two files containing segment and vertex information.

# General keyword arguments
- `dist_threshold::Real`: The distance below which branchpoints are considered to be the same and will be merged.
- `reverse_y::Bool`: Reverse the y-coordinates? Use this to ensure the root system is oriented with the highest y-values at the top (affects results of solving).
- `pₛ::Real`: The quantile of all widths in the graph to use as the width of a single root.
- `nₕ_min::Integer`: The minimum amount of roots possibly present in each edge of the graph.
# Column names
For the file containing vertex/node information
- `node_id_colname`: ID.
- `segment_ids_colname`: IDs of the connected segments.
- `x_colname`: x-coordinate.
- `y_colname`: y-coordinate.
- `lateral_score_colname`: NN-predicted probability of a lateral root dividing at this node.

For the file containing edge/segment information.
- `segment_id_colname`: ID.
- `dist_colname`: Estimated width of segment.
- `primary_score_colname`: NN-predicted probability of the segment containing a primary root.
- `coords_colname`: The y and x coordinates of each pixel of the segment.
"""
function get_supergraph(
        filename_segments::String, filename_vertices::String;
        dist_threshold::Real, reverse_y::Bool,
        node_id_colname = :Node, segment_ids_colname = :Segment_IDs,
        x_colname = :Coord_x, y_colname = :Coord_y, lateral_score_colname = :Lateral_Score,
        segment_id_colname = :Segment_ID, dist_colname = :Mean_Distance,
        primary_score_colname = :Heatmap_Mean, coords_colname = :Coords
    )

    pg = get_pregraph(
        filename_segments, filename_vertices; dist_threshold, reverse_y, node_id_colname,
        segment_ids_colname, x_colname, y_colname, lateral_score_colname, segment_id_colname,
        dist_colname, primary_score_colname, coords_colname
    )
    sg = get_supergraph(pg)

    return sg
end

function get_supergraph(pg::PreGraph)
    Vₕ₀ = [
        HyperVertex(id(mv), HyperEdge.(segments(pg, mv)), x(mv), y(mv), pred_split(mv))
        for mv in getmetavertices(pg)
        if !isspecial(mv)
    ]
    Vₕ₊ = [
        HyperVertex(id(mv), HyperEdge.(segments(pg, mv)), x(mv), y(mv), NaN, [id(mv)])
        for mv in getmetavertices(pg)
        if isspecial(mv)
    ]
    V₀ = [
        [
            SingularVertex(
                v, get_externaledges(v, hv, [Vₕ₀; Vₕ₊]), 
                [InternalEdge(v, roommate(v)), InternalEdge(roommate(v), v)], hv
            )
            for v in vertices(hv)
        ]
        for hv in Vₕ₀
    ] |> x -> reduce(vcat, x)
    V₊ = [
        SingularVertex(id(hv), get_externaledges(hv, [Vₕ₀; Vₕ₊]), InternalEdge{typeof(id(hv))}[], hv)
        for hv in Vₕ₊
    ]

    return SuperGraph(Vₕ₀, Vₕ₊, V₀, V₊)
end

# standard vertices
function get_externaledges(v::T, hv::HyperVertex{T, U}, Vₕ::Vector{HyperVertex{T, U}}) where {T, U}
    externaledges = [
        [
            ExternalEdge(v, v_nb, he)
            for v_nb in vertices(neighbor(hv, he, Vₕ))
            if polarity(v_nb) == polarity(v) || is_augmented(v_nb)
        ]
        for he in edges(hv)
    ] |> x -> reduce(vcat, x)

    return externaledges
end

# augmented vertices
function get_externaledges(hv::HyperVertex{T, U}, Vₕ::Vector{HyperVertex{T, U}}) where {T, U}
    return [
        [
            ExternalEdge(id(hv), v_nb, he)
            for v_nb in vertices(neighbor(hv, he, Vₕ))
        ]
        for he in edges(hv)
    ] |> x -> reduce(vcat, x)
end