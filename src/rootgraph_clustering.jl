"""
    get_subgraphs(rg::RootGraph; augmented_margins::Real = 0.1, pₛ = 0.2, nₕ_min = 1)

Separate all unconnected subgraphs of the graph.

`augmented_margins` controls the position of the augmented vertices on a plot and is purely for aesthetics.
`pₛ` and `nₕ_min` have the same definition as in [`get_rootgraph`](@ref).
"""
function get_subgraphs(rg::RootGraph; augmented_margins::Real = 0.1, pₛ = 0.2, nₕ_min = 1)
    return [
        get_rootgraph(subpregraph(rg, cluster; augmented_margins); pₛ, nₕ_min)
            for cluster in vertex_clusters(rg)
    ]
end

# group the original vertices into connected clusters
function vertex_clusters(rg::RootGraph)
    rvs = copy(V₀(rg))
    clusters = typeof(rvs)[]

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
        push!(clusters, current_cluster)
    end

    return clusters
end

# recreate a cluster of connected vertices as a standalone preliminary graph with fresh consecutive ids
function subpregraph(rg::RootGraph{T, U}, cluster::Vector{RootVertex{T, U}}; augmented_margins) where {T, U}
    # one branchpoint per distinct position in the cluster
    representatives = unique(coords, cluster) # the first vertex at every position
    position2bp = Dict(coords(rv) => convert(T, i) for (i, rv) in enumerate(representatives))
    branchpoints = [
        BranchPoint(position2bp[coords(rv)], Symbol(:bp, position2bp[coords(rv)]), x(rv), y(rv), pred_split(rv))
            for rv in representatives
    ]

    # one pregraph segment per group of parallel edges inside the cluster
    # (a cluster is a connected component, so a segment touching it lies fully inside it)
    cluster_vs = id.(cluster)
    cluster_segments = [seg for seg in segments₀(rg) if first(srcs(seg)) in cluster_vs]
    bp_of(v) = position2bp[coords(getvertex(rg, v))]
    _segments = [
        Segment(convert(T, i), [bp_of(first(srcs(seg))), bp_of(first(dsts(seg)))], width(seg), pred_primary(seg))
            for (i, seg) in enumerate(cluster_segments)
    ]

    augment!(branchpoints, _segments; augmented_margins)

    return get_pregraph(branchpoints, _segments)
end
