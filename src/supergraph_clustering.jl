"""
    get_subgraphs(rg::RootGraph; augmented_margins::Real = 0.1, pₛ = 0.2, nₕ_min = 1)

Separate all unconnected subgraphs of the graph. 

`augmented_margins` controls the position of the augmented vertices on a plot and is purely for aesthetics.
`pₛ` and `nₕ_min` have the same definition as in [`get_supergraph`](@ref).
"""
function get_subgraphs(rg::RootGraph; augmented_margins::Real = 0.1, pₛ = 0.2, nₕ_min = 1)
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
        res = edges.(rv_cluster) |> x -> reduce(vcat, x) |> unique
        all_widths = [width(re) for re in res if !is_augmented(re)]
        single_width = quantile(all_widths, pₛ)
        nₕs = [get_num_hypotheses(rv, single_width; nₕ_min) for rv in rv_cluster]

        id_conversion_dict = Dict([id.(rv_cluster); -3:-1] .=> [eachindex(rv_cluster); -3:-1])

        recreated_V₀ = [recreate_V₀(rv, id_conversion_dict, nₕs, i) for (i, rv) in enumerate(rv_cluster)]
        recreated_Vₕ₊ = [recreate_Vₕ₊(rv, id_conversion_dict, recreated_V₀; augmented_margins) for rv in Vₕ₊(rg)]
        subgraph = RootGraph(
            recreated_V₀,
            recreated_Vₕ₊,
            [
                getsingularvertices(rv, [recreated_V₀; recreated_Vₕ₊]) for rv in recreated_V₀
            ] |> x -> reduce(vcat, x),
            [
                getsingularvertices(rv, [recreated_V₀; recreated_Vₕ₊]) for rv in recreated_Vₕ₊
            ] |> x -> reduce(vcat, x)
        )
        subgraphs[i] = subgraph
    end

    return subgraphs
end

function recreate_V₀(rv::RootVertex{T, U}, id_conversion_dict::Dict, nₕs::Vector{<:Integer}, i::Integer) where {T, U}
    id_new = id_conversion_dict[id(rv)]
    res_new = RootEdge{T, U}[
        RootEdge([id_conversion_dict[v] for v in vertices(re)]..., segment_id(re), width(re), pred_primary(re))
            for re in E(rv) if all(haskey.([id_conversion_dict], vertices(re)))
    ]

    prev_id = sum(nₕs[1:(i - 1)]) # amount of vertices that have been defined in previous hypervertices
    vertices_new = collect(prev_id .+ (1:nₕs[i]))

    rv_new = RootVertex(id_new, res_new, coords(rv)..., pred_split(rv), vertices_new)

    return rv_new
end

function recreate_Vₕ₊(rv::RootVertex{T, U}, id_conversion_dict::Dict, recreated_V₀::Vector{RootVertex{T, U}}; augmented_margins) where {T, U}
    id_new = id_conversion_dict[id(rv)]
    res_new = RootEdge{T, U}[
        RootEdge([id_conversion_dict[v] for v in vertices(re)]..., segment_id(re), width(re), pred_primary(re))
            for re in E(rv) if all(haskey.([id_conversion_dict], vertices(re)))
    ]
    coords_new = get_augmented_coords(id(rv), recreated_V₀; augmented_margins)

    rv_new = RootVertex(id_new, res_new, coords_new..., NaN, vertices(rv))

    return rv_new
end