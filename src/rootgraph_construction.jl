"""
    get_rootgraph(filename_segments::String, filename_vertices::String;
        dist_threshold::Real, reverse_y::Bool, pₛ::Real = 0.2, nₕ_min::Integer = 1, [colnames...]
    )

Create a root graph from two files containing segment and vertex information.

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
function get_rootgraph(
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
    rg = get_rootgraph(pg)

    return rg
end

function get_rootgraph(pg::PreGraph{T, U, V}) where {T, U, V}
    bps = getbranchpoints(pg)

    V₀ = [
        RootVertex(id(bp), [RootEdge(s) for s in segments(pg, bp)], x(bp), y(bp), pred_split(bp))
            for bp in bps if !isspecial(bp)
    ]
    V₊ = [
        RootVertex(id(bp), [RootEdge(s) for s in segments(pg, bp)], x(bp), y(bp), pred_split(bp))
            for bp in bps if isspecial(bp)
    ]

    return RootGraph(V₀, V₊)
end

# create edge from segment
RootEdge(s::Segment) = RootEdge(vertices(s)..., width(s), pred_primary(s))

# get amount of hypotheses corresponding to a branchpoint
# function get_num_hypotheses(pg::PreGraph, bp::BranchPoint; pₛ, nₕ_min)
#     all_widths = [width(s) for s in segments(pg) if !isspecial(s)]
#     single_width = quantile(all_widths, pₛ)
#     nₕ = [nₕ_min + floor(Int64, width(s) / single_width) for s in segments(pg, bp) if !isspecial(s)] |> maximum

#     return nₕ
# end
