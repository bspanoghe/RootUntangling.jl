# # main function

function get_pregraph(edge_data_dict::Dict, vertex_data_dict::Dict; dist_threshold = 3, augmented_margins = 0.3)
    # get metavertices and segments
    metavertices = getmetavertices(vertex_data_dict);
    segments = getsegments(metavertices, vertex_data_dict, edge_data_dict);

    # clean
    differentiate_duplicates!(metavertices, segments)
    remove_single_vertex_segments!(segments)
    cluster_vertices!(segments, metavertices; dist_threshold);
    remove_unconnected_vertices!(metavertices, segments)
    remake_ids!(metavertices, segments)

    # augmentation
    augment!(metavertices, segments; augmented_margins)
    
    # turn into preliminary graph
    pg = get_pregraph(metavertices, segments);
    
    return pg
end

# # get metavertices & segments from file
# ## metavertices
function getmetavertices(vertex_data_dict::Dict)
    metavertices = MetaVertex{Int64, Float64}[
        MetaVertex(i, Symbol(vertex_entry[1]), [float(vertex_entry[2][col]) for col in [:x, :y, :pred_split]]...)
        for (i, vertex_entry) in enumerate(vertex_data_dict)
    ]
    return metavertices
end

# ## segments
function getsegments(metavertices::Vector{MetaVertex{T, U}}, vertex_data_dict, edge_data_dict) where {T, U}
    segment_connections = Dict{T, Vector{T}}() # edge id => vertices
    for metavertex in metavertices
        vertex_data = vertex_data_dict[metavertex.name]

        for edge_id in vertex_data[:segment_ids]
            vertices = get(segment_connections, edge_id, T[])
            segment_connections[edge_id] = [vertices; metavertex.id] 
        end
    end

    segments = [
        Segment(id, vertices, [edge_data_dict[id][col] for col in [:width, :pred_primary, :xs, :ys]]...)
        for (id, vertices) in segment_connections
    ]

    return segments
end

# # cleaning

# ## break segments with identical vertices in two
# adds new metavertices and segments to differentiate them
# like so                                                
#      /---\            /-o-\                              
#   --o    o--   =>  --o    o--                                
#     \---/            \-o-/                              
function differentiate_duplicates!(metavertices::Vector{MetaVertex{T, U}},
    segments::Vector{Segment{T, V}}) where {T, U, V}

    # find "duplicate" segments
    duplicates = duplicate_elements(s -> sort(vertices(s)), segments)
    isempty(duplicates) && return nothing

    # make new metavertices at halfpoints of duplicates
    v_max = id.(metavertices) |> maximum
    seg_id_max = id.(segments) |> maximum

    mv_new = MetaVertex{T, U}[
        MetaVertex(
            v_max+i,
            Symbol("hp_$(i)"),
            xs(duplicates[i])[end ÷ 2] |> x -> convert(U, x),
            ys(duplicates[i])[end ÷ 2] |> x -> convert(U, x),
            zero(U)
        )
        for i in eachindex(duplicates)
    ]

    # make new segments between startpoint and new midpoint
    seg_new = [
        Segment{T, V}(
            seg_id_max + i,
            [vertices(duplicates[i])[1], id(mv_new[i])],
            [f(duplicates[i]) for f in [width, pred_primary, xs, ys]]...
        )
        for i in eachindex(duplicates)
    ]

    # replace startpoint with new halfpoint
    for i in eachindex(duplicates)
        vertices(duplicates[i])[1] = id(mv_new[i])
    end

    # add everything to variables
    append!(metavertices, mv_new)
    append!(segments, seg_new)
    
    return nothing
end

function duplicate_elements(v::Vector)
    seen = Dict{eltype(v), Ref{Int}}() # ty julia discourse user oxinabox, and for making me read about Refs
    for x in v
        get!(() -> 0, seen, x)[] += 1
    end
    
    return [x for x in v if seen[x][] > 1]
end

function duplicate_elements(f::Function, v::Vector)
    v_id = f.(v)

    seen = Dict{eltype(v_id), Ref{Int}}()
    for x in v_id
        get!(() -> 0, seen, x)[] += 1
    end
    
    return [v[i] for i in eachindex(v) if seen[v_id[i]][] > 1]
end

# ## segment cleaning
function remove_single_vertex_segments!(segments)
    bad_segments_idxs = [length(vertices(s)) == 1 for s in segments] |> findall
    !isempty(bad_segments_idxs) && @info "$(length(bad_segments_idxs)) segments found connected to only a single vertex. This is a masking artefact and may usually be safely ignored."
    deleteat!(segments, bad_segments_idxs)

    return nothing
end

# ## vertex clustering
function cluster_vertices!(segments, metavertices; dist_threshold)
    vertex_clusters = get_vertex_clusters(segments, metavertices, dist_threshold)

    if isempty(vertex_clusters)
        return nothing
    end
    check_clusters(vertex_clusters)
    @info "The following clusters were found: $([name(mv) for mv in metavertices if id(mv) in reduce(vcat, vertex_clusters)])"

    merged_metavertices = [
        make_merged_metavertex(vertex_cluster, metavertices, length(metavertices) + i)
        for (i, vertex_cluster) in enumerate(vertex_clusters)
    ]
    append!(metavertices, merged_metavertices)
    merge_clusters!(segments, vertex_clusters, merged_metavertices)

    remove_bad_eps!(segments, metavertices)
    remove_clusternodes!(metavertices, vertex_clusters) # remove clustered vertices after checking segments so we can use their names in the log

    return nothing
end

# ### collect vertices that should be clustered
function get_vertex_clusters(segments::Vector{Segment{T, U}}, metavertices::Vector{<:MetaVertex{T}}, dist_threshold) where {T, U}
    ov_vertices = [vertices(s) for s in segments if length(vertices(s)) != 2] |>
        x -> reduce(vcat, x, init = T[]) |> unique # overconnected vertex ids

    vertex_clusters = Vector{T}[]
    for ov_vertex in ov_vertices
        cluster = T[]

        connected_segments = [s for s in segments if ov_vertex in vertices(s)]
        connected_vertices = vertices.(connected_segments)

        for con_vertex in unique(reduce(vcat, connected_vertices)) 
            if distance(metavertices[ov_vertex], metavertices[con_vertex]) <= dist_threshold
                # connect all vertices that are connected by at least one segment and are really close to the main vertex
                push!(cluster, con_vertex)
            end
        end

        length(cluster) <= 1 && continue # skip rest if cluster has one or no vertices

        same_cluster_idxs = findall.([in.(v, vertex_clusters) for v in cluster]) |> 
            x -> reduce(vcat, x) |> unique
        if !isempty(same_cluster_idxs) # other cluster exists that shares vertices from current cluster
            same_cluster = vertex_clusters[only(same_cluster_idxs)] # should be only maximum 1 cluster with shared vertices
            ov_vertex in same_cluster || push!(same_cluster, ov_vertex) 
                # if there's clusters containing any other vertices from the current cluster, treat it as the same cluster
        else
            push!(vertex_clusters, cluster)
        end
    end

    return vertex_clusters
end

# ### check if vertex clusters make sense
function check_clusters(vertex_clusters)
    cluster_vertices = reduce(vcat, vertex_clusters)
    @assert length(cluster_vertices) == length(unique(cluster_vertices))
end

# ### add a merged metavertex to the metavertices, based on a clusters of vertices
function make_merged_metavertex(vertex_cluster::Vector{T}, metavertices::Vector{<:MetaVertex{T}}, cluster_id) where {T}
    cluster_id = convert(T, cluster_id)
    vertex_names = sort([name(metavertices[vertex]) for vertex in vertex_cluster])
    cluster_name = [string(nn) * (nn == vertex_names[end] ? "" : "_") for nn in vertex_names] |> x -> *(x...) |> Symbol
    cluster_stats = [my_mean([stat(metavertices[v]) for v in vertex_cluster]) for stat in [x, y, pred_split]]
    
    merged_metavertex = MetaVertex(cluster_id, cluster_name, cluster_stats...)

    return merged_metavertex
end
my_mean(xs::AbstractArray{T}) where {T} = convert(T, sum(xs) / length(xs)) # type stable mean
my_mean(xs::AbstractArray{T}) where {T <: Integer} = round(T, sum(xs) / length(xs))

# ### replace vertices from vertex clusters with merged vertices
function merge_clusters!(segments, vertex_clusters, merged_metavertices)
    for s in segments
        for (cluster, merged_metavertex) in zip(vertex_clusters, merged_metavertices)
            in_cluster_idxs = [vertex in cluster for vertex in vertices(s)] |> findall
            if !isempty(in_cluster_idxs)
                deleteat!(s.vertices, in_cluster_idxs)
                append!(s.vertices, merged_metavertex.id)
            end
        end
    end

    return nothing
end

# ### faulty endpoint removal: for segments connected with 2 bps, remove any potential connected endpoints
function remove_bad_eps!(segments, metavertices)
    overconnected_segments = [s for s in segments if length(vertices(s)) > 2]
    for overconnected_segment in overconnected_segments
        segment_metavertices = metavertices[vertices(overconnected_segment)]
        are_endpoints = [name(v)[1] == 'e' for v in segment_metavertices]
        are_branchpoints = .!are_endpoints

        if sum(are_branchpoints) >= 2
            deleteat!(vertices(overconnected_segment), findall(are_endpoints))
        end
    end
    
    return nothing
end

# ### remove the original vertices that have been replaced by a clustered vertex
function remove_clusternodes!(metavertices, vertex_clusters)
    deleteat!(metavertices, sort(unique(reduce(vcat, vertex_clusters))))

    return nothing
end

# ## vertex cleaning
function remove_unconnected_vertices!(metavertices, segments)
    segment_vertices = reduce(vcat, vertices.(segments))
    unconnected_vertex_idxs = findall(mv -> !any(id(mv) in segment_vertices), metavertices)
    deleteat!(metavertices, unconnected_vertex_idxs)

    return nothing
end

# ## id restructuring 
# the ids of the non-augmented vertices should equal `1:length(vertices)`
function remake_ids!(metavertices, segments)
    id_conversion_dict = Dict([
        id(mv) => i
        for (i, mv) in enumerate(metavertices)
    ])

    for (i, mv) in enumerate(metavertices)
        metavertices[i].id = id_conversion_dict[id(mv)]
    end

    for (i, seg) in enumerate(segments)
        segments[i].vertices = [id_conversion_dict[v] for v in vertices(seg)]
    end

    return nothing
end

# # augmentation
function augment!(metavertices, segments; augmented_margins)
    pg₀ = get_pregraph(metavertices, segments)

    v_a = MetaVertex(-1, :appearance, get_augmented_coords(-1, metavertices; augmented_margins)...)
    v_d = MetaVertex(-2, :disappearance, get_augmented_coords(-2, metavertices; augmented_margins)...)
    v_s = MetaVertex(-3, :splitting, get_augmented_coords(-3, metavertices; augmented_margins)...)

    appearance_segments = [Segment(-1, [-1, i]) for i in id.(metavertices)]
    disappearance_segments = [Segment(-2, [-2, i]) for i in id.(metavertices)]
    splitting_segments = [Segment(-3, [-3, i]) for i in [v for v in id.(metavertices) if v in inner_vertices(pg₀)]]

    append!(metavertices, [v_a, v_d, v_s])
    append!(segments, [appearance_segments; disappearance_segments; splitting_segments])

    return nothing
end

function get_augmented_coords(id, vertices; augmented_margins)
    @assert id ∈ -1:-1:-3 "Augmented vertex positions were defined for 3 augmented vertices only"
    ymin, ymax = extrema(y.(vertices))
    xmin, xmax = extrema(x.(vertices))
    x_vertex = xmin - augmented_margins*(xmax-xmin)
    y_vertex = ymax - (abs(id)-1)/2 * (ymax-ymin)

    return x_vertex, y_vertex
end

# # turn into preliminary graph
function get_pregraph(metavertices::Vector{MetaVertex{T, U}}, segments::Vector{Segment{T, V}}) where {T, U, V}
    _vertices = id.(metavertices)
    _segments = [s for s in segments if length(vertices(s)) == 2]
    length(segments) != length(_segments) && @info "Segments were found connecting more than 2 vertices. These have been removed."
    _metavertexdict = Dict(Pair.(_vertices, metavertices))

    return PreGraph(_vertices, _segments, _metavertexdict)
end

# # method for doing all steps from file reading
function get_pregraph(filename_segments::String, filename_vertices::String;
        dist_threshold::Real, reverse_y::Bool, node_id_colname = :Node, segment_ids_colname = :Segment_IDs,
        x_colname = :Coord_x, y_colname = :Coord_y, lateral_score_colname = :Lateral_Score,
        segment_id_colname = :Segment_ID, dist_colname = :Mean_Distance,
        primary_score_colname = :Heatmap_Mean, coords_colname = :Coords
    )

    vertex_data = read_data(filename_vertices, node_id_colname)
    ymax = [vertex_datum[y_colname] for vertex_datum in values(vertex_data)] |> maximum
    y_transform = reverse_y ? y -> ymax .- y : identity
    vertex_data_dict = get_vertex_info(vertex_data; segment_ids_colname, x_colname,
        y_colname, lateral_score_colname, y_transform
    )

    edge_data = read_data(filename_segments, segment_id_colname)
    edge_data_dict = get_edge_info(edge_data; dist_colname,
        primary_score_colname, coords_colname, y_transform
    )

    pg = get_pregraph(edge_data_dict, vertex_data_dict; dist_threshold)

    return pg
end