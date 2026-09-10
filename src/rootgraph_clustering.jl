"""
    get_subgraphs(rg::RootGraph; augmented_margins::Real = 0.1)

Separate all unconnected subgraphs of the graph. 

`augmented_margins` controls the position of the augmented vertices on a plot and is purely for aesthetics.
"""
function get_subgraphs(rg::RootGraph; augmented_margins::Real = 0.1)
    rvs = copy(V₀(rg))
    rv_clusters = typeof(rvs)[]

    # cluster connected hypervertices
    while !isempty(rvs)
        current_cluster = rvs[1:1]
        deleteat!(rvs, 1)
        i = 0
        while length(current_cluster) > i
            i += 1
            for nb in neighbors(rg, current_cluster[i])
                nb_idx = findfirst(x -> x == nb, rvs)
                if !isnothing(nb_idx)
                    push!(current_cluster, nb)
                    deleteat!(rvs, nb_idx)
                end
            end
        end
        push!(rv_clusters, current_cluster)
    end

    # turn rv clusters into RootGraphs
    subgraphs = Vector{typeof(rg)}(undef, length(rv_clusters))
    for (i, rv_cluster) in enumerate(rv_clusters)
        id_conversion_dict = Dict([id.(rv_cluster); -3:-1] .=> [eachindex(rv_cluster); -3:-1])

        recreated_V₀ = [recreate_V₀(rv, id_conversion_dict) for rv in rv_cluster]
        recreated_V₊ = [recreate_V₊(rv, id_conversion_dict, recreated_V₀; augmented_margins) for rv in V₊(rg)]
        subgraph = RootGraph(
            recreated_V₀,
            recreated_V₊,
        )
        subgraphs[i] = subgraph
    end

    return subgraphs
end

function recreate_V₀(rv::RootVertex{T, U}, id_conversion_dict::Dict) where {T, U}
    id_new = id_conversion_dict[id(rv)]
    res_new = RootEdge{T, U}[
        RootEdge([id_conversion_dict[v] for v in vertices(re)]..., segment_id(re), width(re), pred_primary(re))
        for re in E(rv) if all(haskey.([id_conversion_dict], vertices(re)))
    ]

    rv_new = RootVertex(id_new, res_new, coords(rv)..., pred_split(rv))

    return rv_new
end

function recreate_V₊(rv::RootVertex{T, U}, id_conversion_dict::Dict, recreated_V₀::Vector{RootVertex{T, U}}; augmented_margins) where {T, U}
    id_new = id_conversion_dict[id(rv)]
    res_new = RootEdge{T, U}[
        RootEdge([id_conversion_dict[v] for v in vertices(re)]..., segment_id(re), width(re), pred_primary(re))
        for re in E(rv) if all(haskey.([id_conversion_dict], vertices(re)))
    ]
    coords_new = get_augmented_coords(id(rv), recreated_V₀; augmented_margins)

    rv_new = RootVertex(id_new, res_new, coords_new..., NaN)

    return rv_new
end