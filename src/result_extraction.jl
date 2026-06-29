# # classification of graph edges
# one graph

# function get_se_classification_dict(sg::SuperGraph, model::JuMP.Model)
#     eₚ = [var for var in all_variables(model) if !isnothing(match(r"^eₚ\[\d+\]$", JuMP.name(var)))]
#     eₗ = [var for var in all_variables(model) if !isnothing(match(r"^eₗ\[\d+\]$", JuMP.name(var)))]

#     se_classification_dict = [E(sg)[i] => round(Bool, value(eₚ[i])) + round(Bool, value(eₗ[i]))*im for i in eachindex(E(sg))] |> Dict

#     return se_classification_dict
# end

function get_se_classification_dict(sg::SuperGraph, model::JuMP.Model)
    e = [var for var in all_variables(model) if !isnothing(match(r"^e\[\d+\]$", JuMP.name(var)))]
    eₚ = [var for var in all_variables(model) if !isnothing(match(r"^eₚ\[\d+\]$", JuMP.name(var)))]

    se_classification_dict = [E(sg)[i] => round(Bool, value(eₚ[i])) + round(Bool, value(e[i]) - value(eₚ[i]))*im for i in eachindex(E(sg))] |> Dict

    return se_classification_dict
end

function get_he_classification_dict(sg::SuperGraph, model::JuMP.Model)
    se_classification_dict = get_se_classification_dict(sg, model)
    he_classification_dict = [he => sum([se_classification_dict[e] for e in E(sg, he)]) for he in Eₕ(sg)] |> Dict

    return he_classification_dict
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
    vertex_ids::Vector{T}
    xs::Vector{U}
    ys::Vector{U}
    function Root(is_primary, vertex_ids, xs, ys)
        @assert length(vertex_ids) == length(xs) == length(ys)
        return new{eltype(vertex_ids), eltype(xs)}(is_primary, vertex_ids, xs, ys)
    end
end
is_primary(r::Root) = r.is_primary
vertex_ids(r::Root) = r.vertex_ids
xs(r::Root) = r.xs
ys(r::Root) = r.ys

function Root(
        se::SingularEdge{T, U}, se_classification_dict::Dict{SingularEdge{T, U}, Complex{Int64}}, sg::SuperGraph{T, U}
    ) where {T, U}
    @assert haskey(se_classification_dict, se)
    Root(
        real(se_classification_dict[se]) == 1,
        [vertices(se)...],
        xs(sg, se),
        ys(sg, se)
    )
end

function get_roots(sg::SuperGraph, model::JuMP.Model)
    se_classification_dict = get_se_classification_dict(sg, model)
    active_edges = [e for e in E₀(sg) if abs(se_classification_dict[e]) > 0]
    
    roots = Root[]
    while !isempty(active_edges)
        current_root = Root(active_edges[1], se_classification_dict, sg)
        deleteat!(active_edges, 1)
        edge_idx = findfirst(x -> are_connected(current_root, x), active_edges)
        while !isnothing(edge_idx)
            push!(current_root, active_edges[edge_idx], sg)
            deleteat!(active_edges, edge_idx)
            edge_idx = findfirst(x -> are_connected(current_root, x), active_edges)
        end
        push!(roots, current_root)
    end

    return roots
end

function Base.push!(r::Root{T, U}, se::SingularEdge{T, U}, sg::SuperGraph{T, U}) where {T, U}
    new_vertex_idx = findfirst(x -> !(x in vertex_ids(r)[[1, end]]), vertices(se))
    if isnothing(new_vertex_idx)
        @warn "Loop found in root"
        new_vertex = vertex_ids(r)[1]
    else
        new_vertex = vertices(se)[new_vertex_idx]
    end
    new_sv = getsingularvertex(sg, new_vertex)

    is_upstream = vertex_ids(r)[1] in vertices(se)
    if is_upstream
        pushfirst!(vertex_ids(r), new_vertex)
        pushfirst!(xs(r), x(new_sv))
        pushfirst!(ys(r), y(new_sv))
    else
        push!(vertex_ids(r), new_vertex)
        push!(xs(r), x(new_sv))
        push!(ys(r), y(new_sv))
    end

    return nothing
end
are_connected(r::Root{T, U}, se::SingularEdge{T, U}) where {T, U} = (
    !isempty(intersect(vertex_ids(r)[[1, end]], vertices(se)))
)