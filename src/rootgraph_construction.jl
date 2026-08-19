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
        dist_threshold::Real, reverse_y::Bool, pₛ::Real = 0.2, nₕ_min::Integer = 1,
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
    rg = get_rootgraph(pg; pₛ, nₕ_min)

    return rg
end

function get_rootgraph(pg::PreGraph{T, U, V}; pₛ, nₕ_min) where {T, U, V}
    bps = getbranchpoints(pg)

    # assign each branchpoint its vertex ids:
    # one hypothetical root vertex per possible root passing through an original branchpoint,
    # the branchpoint id itself for augmented branchpoints
    bp2vs = Dict{T, Vector{T}}()
    v_max = 0
    for bp in bps
        nₕ = isspecial(bp) ? 0 : get_num_hypotheses(pg, bp; pₛ, nₕ_min)
        bp2vs[id(bp)] = isspecial(bp) ? [id(bp)] : collect(v_max .+ (1:nₕ))
        v_max += nₕ
    end

    # connect each vertex of a branchpoint to all vertices of its neighboring branchpoints
    rootedges(bp, v) = [
        RootEdge(v, w, width(s), pred_primary(s))
            for s in segments(pg, bp)
            for w in bp2vs[other_vertex(s, id(bp))]
    ]

    V₀ = [
        RootVertex(v, rootedges(bp, v), x(bp), y(bp), pred_split(bp))
            for bp in bps if !isspecial(bp)
            for v in bp2vs[id(bp)]
    ]
    V₊ = [
        RootVertex(id(bp), rootedges(bp, id(bp)), x(bp), y(bp), NaN)
            for bp in bps if isspecial(bp)
    ]

    return RootGraph(V₀, V₊)
end

# the vertex at the other end of a segment
other_vertex(s::Segment{T, U}, v::T) where {T, U} = vertices(s)[1] == v ? vertices(s)[2] : vertices(s)[1]

# get amount of hypotheses corresponding to a branchpoint
function get_num_hypotheses(pg::PreGraph, bp::BranchPoint; pₛ, nₕ_min)
    all_widths = [width(s) for s in segments(pg) if !isspecial(s)]
    single_width = quantile(all_widths, pₛ)
    nₕ = [nₕ_min + floor(Int64, width(s) / single_width) for s in segments(pg, bp) if !isspecial(s)] |> maximum

    return nₕ
end
