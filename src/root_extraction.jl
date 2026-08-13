# # classification of graph edges/vertices
# one graph

function get_se_classification_dict(sg::SuperGraph, model::JuMP.Model)
    ea = [var for var in all_variables(model) if !isnothing(match(r"^ea\[\d+\]$", JuMP.name(var)))]
    ep = [var for var in all_variables(model) if !isnothing(match(r"^ep\[\d+\]$", JuMP.name(var)))]

    se_classification_dict = [E(sg)[i] => round(Bool, value(ep[i])) + round(Bool, value(ea[i]) - value(ep[i]))*im for i in eachindex(E(sg))] |> Dict

    return se_classification_dict
end

function get_he_classification_dict(sg::SuperGraph, model::JuMP.Model)
    se_classification_dict = get_se_classification_dict(sg, model)
    he_classification_dict = [he => sum([se_classification_dict[e] for e in E(sg, he)]) for he in Eₕ(sg)] |> Dict

    return he_classification_dict
end

function get_sv_classification_dict(sg::SuperGraph, model::JuMP.Model)
    va = [var for var in all_variables(model) if !isnothing(match(r"^va\[\d+\]$", JuMP.name(var)))]
    vp = [var for var in all_variables(model) if !isnothing(match(r"^vp\[\d+\]$", JuMP.name(var)))]

    sv_classification_dict = [V₀(sg)[i] => round(Bool, value(vp[i])) + round(Bool, value(va[i]) - value(vp[i]))*im for i in eachindex(V₀(sg))] |> Dict

    return sv_classification_dict
end

function get_hv_classification_dict(sg::SuperGraph, model::JuMP.Model)
    sv_classification_dict = get_sv_classification_dict(sg, model)
    hv_classification_dict = [hv => sum([sv_classification_dict[v] for v in V(sg, hv)]) for hv in Vₕ₀(sg)] |> Dict

    return hv_classification_dict
end

# multiple graphs

function get_se_classification_dict(sgs::Vector{SuperGraph{T, U}}, models::Vector{<:JuMP.Model}; kwargs...) where {T, U}
    classification_dict = Dict{Int64, Dict{SingularEdge{T, U}, Complex{Int64}}}()

    for (i, (sg, model)) in enumerate(zip(sgs, models))
        classification_dict[i] = get_se_classification_dict(sg, model)
    end
    
    return classification_dict
end

function get_he_classification_dict(sgs::Vector{SuperGraph{T, U}}, models::Vector{<:JuMP.Model}; kwargs...) where {T, U}
    classification_dict = Dict{Int64, Dict{HyperEdge{T, U}, Complex{Int64}}}()

    for (i, (sg, model)) in enumerate(zip(sgs, models))
        classification_dict[i] = get_he_classification_dict(sg, model)
    end
    
    return classification_dict
end

# # extract roots

struct Root{T, U}
    is_primary::Bool
    V::Vector{SingularVertex{T, U}}
end
is_primary(r::Root) = r.is_primary
V(r::Root) = r.V

vertices(r::Root) = id.(V(r))
xs(r::Root) = x.(V(r))
ys(r::Root) = y.(V(r))

Base.length(r::Root) = length(r.V)

function Root(
        se::SingularEdge{T, U}, se_classification_dict::Dict{SingularEdge{T, U}, Complex{Int64}},
        sg::SuperGraph{T, U}
    ) where {T, U}
    @assert haskey(se_classification_dict, se)
    Root(
        real(se_classification_dict[se]) == 1,
        [getsingularvertex(sg, v) for v in vertices(se)],
    )
end

function get_roots(sg::SuperGraph, model::JuMP.Model)
    se_classification_dict = get_se_classification_dict(sg, model)
    active_edges = [e for e in E₀(sg) if abs(se_classification_dict[e]) > 0]
    
    roots = Root[]
    while !isempty(active_edges)
        current_root = Root(active_edges[1], se_classification_dict, sg)
        deleteat!(active_edges, 1)

        edge_idx = findfirst(e -> are_connected(current_root, e), active_edges)
        while !isnothing(edge_idx)
            grow!(current_root, active_edges[edge_idx], sg)
            deleteat!(active_edges, edge_idx)
            edge_idx = findfirst(e -> are_connected(current_root, e), active_edges)
        end

        push!(roots, current_root)
    end

    sort!(roots, by = r -> is_primary(r), rev = true)

    return roots
end
are_connected(r::Root{T, U}, se::SingularEdge{T, U}) where {T, U} = (
    !isempty(intersect(vertices(r)[[1, end]], vertices(se)))
)

function grow!(r::Root{T, U}, se::SingularEdge{T, U}, sg::SuperGraph{T, U}) where {T, U}
    new_vertex_idx = findfirst(x -> !(x in vertices(r)[[1, end]]), vertices(se))
    if isnothing(new_vertex_idx)
        @warn "Loop found in root"
        new_vertex = vertices(r)[1]
    else
        new_vertex = vertices(se)[new_vertex_idx]
    end

    is_upstream = vertices(r)[1] in vertices(se)
    if is_upstream
        pushfirst!(V(r), getsingularvertex(sg, new_vertex))
    else
        push!(V(r), getsingularvertex(sg, new_vertex))
    end

    return nothing
end
