"""
    get_subgraphs(sg::SuperGraph; augmented_margins::Real = 0.1, pₛ = 0.2, nₕ_min = 1)

Separate all unconnected subgraphs of the graph. 

`augmented_margins` controls the position of the augmented vertices on a plot and is purely for aesthetics.
`pₛ` and `nₕ_min` have the same definition as in [`get_supergraph`](@ref).
"""
function get_subgraphs(sg::SuperGraph; augmented_margins::Real = 0.1, pₛ = 0.2, nₕ_min = 1)
    hvs = copy(Vₕ₀(sg))
    hv_clusters = typeof(hvs)[]

    # cluster connected hypervertices
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

    # turn hv clusters into SuperGraphs
    subgraphs = Vector{typeof(sg)}(undef, length(hv_clusters))
    for (i, hv_cluster) in enumerate(hv_clusters)
        id_conversion_dict = Dict([id.(hv_cluster); -3:-1] .=> [eachindex(hv_cluster); -3:-1])
        nₕs = [get_num_hypotheses(hv_cluster, hv; pₛ, nₕ_min) for hv in hv_cluster]

        recreated_Vₕ₀ = [recreate_Vₕ₀(hv, id_conversion_dict, nₕs, i) for (i, hv) in enumerate(hv_cluster)]
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

function recreate_Vₕ₀(hv::HyperVertex{T, U}, id_conversion_dict::Dict, nₕs::Vector{<:Integer}, i::Integer) where {T, U}
    id_new = id_conversion_dict[id(hv)]
    hes_new = HyperEdge{T, U}[
        HyperEdge([id_conversion_dict[v] for v in vertices(he)]..., width(he), pred_primary(he))
            for he in E(hv) if all(haskey.([id_conversion_dict], vertices(he)))
    ]

    prev_id = sum(nₕs[1:(i - 1)]) # amount of vertices that have been defined in previous hypervertices
    vertices_new = collect(prev_id .+ (1:nₕs[i]))

    hv_new = HyperVertex(id_new, hes_new, coords(hv)..., pred_split(hv), vertices_new)

    return hv_new
end

function recreate_Vₕ₊(hv::HyperVertex{T, U}, id_conversion_dict::Dict, recreated_Vₕ₀::Vector{HyperVertex{T, U}}; augmented_margins) where {T, U}
    id_new = id_conversion_dict[id(hv)]
    hes_new = HyperEdge{T, U}[
        HyperEdge([id_conversion_dict[v] for v in vertices(he)]..., width(he), pred_primary(he))
            for he in E(hv) if all(haskey.([id_conversion_dict], vertices(he)))
    ]
    coords_new = get_augmented_coords(id(hv), recreated_Vₕ₀; augmented_margins)

    hv_new = HyperVertex(id_new, hes_new, coords_new..., NaN, vertices(hv))

    return hv_new
end

function get_num_hypotheses(hv_cluster::Vector{<:HyperVertex}, hv::HyperVertex; pₛ, nₕ_min)
    hes = edges.(hv_cluster) |> x -> reduce(vcat, x) |> unique
    all_widths = [width(he) for he in hes if !is_augmented(he)]
    single_width = quantile(all_widths, pₛ)
    nₕ = [nₕ_min + floor(Int64, width(he)/single_width) for he in edges(hv) if !is_augmented(he)] |> maximum

    return nₕ
end
