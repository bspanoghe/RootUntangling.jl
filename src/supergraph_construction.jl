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
    sg = get_supergraph(pg; pₛ, nₕ_min)

    return sg
end

function get_supergraph(pg::PreGraph; pₛ, nₕ_min)
    nₕs = [get_num_hypotheses(pg, mv; pₛ, nₕ_min) for mv in getmetavertices(pg) if !isspecial(mv)]

    Vₕ₀ = [
        get_hypervertex(pg, mv, nₕs, i)
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
function get_num_hypotheses(pg::PreGraph, mv::MetaVertex; pₛ, nₕ_min)
    all_widths = [width(s) for s in segments(pg) if !isspecial(s)]
    single_width = quantile(all_widths, pₛ)
    nₕ = [nₕ_min + floor(Int64, width(s)/single_width) for s in segments(pg, mv) if !isspecial(s)] |> maximum

    return nₕ
end

# instantiate a hypervertex based on a metavertex
function get_hypervertex(pg::PreGraph{T, U, V}, mv::MetaVertex{T, U}, nₕs::Vector{<:Integer}, i::Integer) where {T, U, V}
    prev_id = sum(nₕs[1:(i - 1)]) # amount of vertices that have been defined in previous hypervertices
    vertices = collect(prev_id .+ (1:nₕs[i]))

    return HyperVertex(id(mv), HyperEdge.(segments(pg, mv)), x(mv), y(mv), pred_split(mv), vertices)
end

HyperEdge(s::Segment) = HyperEdge(vertices(s)..., width(s), pred_primary(s))
SingularEdge(s::Segment) = SingularEdge(vertices(s)..., HyperEdge(s))

# get singular vertices corresponding with a given hypervertex (containing the correct edges)
function getsingularvertices(hv::HyperVertex{T, U}, Vₕ::Vector{HyperVertex{T, U}}) where {T, U}
    return [
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
