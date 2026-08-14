# # Type
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
curve_length(r::Root) = sqrt.(diff(xs(r)).^2 + diff(ys(r)).^2) |> sum
distance(r::Root) = sqrt( (ys(r)[end] - ys(r)[1])^2 +  (xs(r)[end] - xs(r)[1])^2 )
tortuosity(r::Root) = curve_length(r) / distance(r)
tortuosity(rs::Vector{<:Root}) = sum(tortuosity.(rs)) / length(rs)

# # Construction
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