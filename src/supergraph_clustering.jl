function get_subgraphs(sg::SuperGraph; augmented_margins::Real = 0.1)
    hvs = copy(Vₕ₀(sg))
    hv_clusters = typeof(hvs)[]

    while !isempty(hvs)
        current_cluster = hvs[1:1]
        deleteat!(hvs, 1)
        i = 0
        while length(current_cluster) > i
            i += 1
            for nb in neighbors(sg, current_cluster[i])
                nb_idx = findfirst(x -> x == nb, hvs)
                if !isnothing(nb_idx)
                    push!(current_cluster, nb)
                    deleteat!(hvs, nb_idx)
                end
            end
        end
        push!(hv_clusters, current_cluster)
    end

    subgraphs = Vector{typeof(sg)}(undef, length(hv_clusters))
    for (i, hv_cluster) in enumerate(hv_clusters)
        id_conversion_dict = Dict([id.(hv_cluster); -3:-1] .=> [eachindex(hv_cluster); -3:-1])

        recreated_Vₕ₀ = [recreate_Vₕ₀(hv, id_conversion_dict) for hv in hv_cluster]
        recreated_Vₕ₊ = [recreate_Vₕ₊(hv, id_conversion_dict, recreated_Vₕ₀; augmented_margins) for hv in Vₕ₊(sg)]
        subgraph = SuperGraph(
            recreated_Vₕ₀,
            recreated_Vₕ₊, 
            [
                getsingularvertices(hv, [recreated_Vₕ₀; recreated_Vₕ₊]) for hv in recreated_Vₕ₀
            ] |> x -> reduce(vcat, x), 
            [
                getsingularvertices(hv, [recreated_Vₕ₀; recreated_Vₕ₊]) for hv in recreated_Vₕ₊
            ] |> x -> reduce(vcat, x)
        )
        subgraphs[i] = subgraph
    end

    return subgraphs
end

function recreate_Vₕ₀(hv::HyperVertex{T, U}, id_conversion_dict::Dict) where {T, U}
    id_new = id_conversion_dict[id(hv)]
    hes_new = HyperEdge{T, U}[
            HyperEdge([id_conversion_dict[v] for v in vertices(he)]..., width(he), pred_primary(he))
            for he in E(hv) if all(haskey.([id_conversion_dict], vertices(he)))
    ]

    hv_new = HyperVertex(id_new, hes_new, coords(hv)..., pred_split(hv), num_hypotheses(hv))

    return hv_new
end

function recreate_Vₕ₊(hv::HyperVertex{T, U}, id_conversion_dict::Dict, recreated_Vₕ₀::Vector{HyperVertex{T, U}}; augmented_margins) where {T, U}
    id_new = id_conversion_dict[id(hv)]
    hes_new = HyperEdge{T, U}[
            HyperEdge([id_conversion_dict[v] for v in vertices(he)]..., width(he), pred_primary(he))
            for he in E(hv) if all(haskey.([id_conversion_dict], vertices(he)))
    ]
    coords_new = get_augmented_coords(id(hv), recreated_Vₕ₀; augmented_margins)

    hv_new = HyperVertex(id_new, hes_new, coords_new..., NaN, num_hypotheses(hv))

    return hv_new
end