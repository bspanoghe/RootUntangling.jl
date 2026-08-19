"""
    get_roots(sg::SuperGraph, model::JuMP.Model)

Extract the roots from a graph `sg` and its solution contained in `model`.

See also [`Root`](@ref).
"""
function get_roots(sg::SuperGraph, model::JuMP.Model)
    se_classification_dict = get_se_classification_dict(sg, model)
    polarity_classification_dict = get_polarity_classification_dict(sg, model)
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
        correct_polarity!(sg, polarity_classification_dict, current_root)
        push!(roots, current_root)
    end

    if length(filter(r -> is_primary(r), roots)) == 1
        sort_root_system!(roots)
        return roots
    else
        return separate_root_systems(sg, se_classification_dict, roots)
    end
end

are_connected(r::Root{T, U}, se::SingularEdge{T, U}) where {T, U} = (
    !isempty(intersect(vertices(r)[[1, end]], vertices(se)))
)

function grow!(r::Root{T, U}, se::SingularEdge{T, U}, sg::SuperGraph{T, U}) where {T, U}
    new_vertex_idx = findfirst(x -> !(x in vertices(r)[[1, end]]), vertices(se))
    if isnothing(new_vertex_idx)
        @warn "Loop found in root"
        return nothing
        # new_vertex = vertices(r)[1]
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

function correct_polarity!(sg::SuperGraph, polarity_classification_dict::Dict, r::Root)
    se = E₀(sg)[findfirst(e -> issetequal(vertices(e), vertices(r)[1:2]), E₀(sg))]
    polarity_match = polarity(V(r)[1:2]...) == polarity_classification_dict[se]
    if !polarity_match
        reverse!(V(r))
    end

    return nothing
end

# divide roots into separate root systems
function separate_root_systems(sg::SuperGraph, se_classification_dict::Dict, rs::Vector{Root})
    primary_roots = filter(is_primary, rs)
    lateral_roots = filter(!is_primary, rs)
    root_systems = [[pr] for pr in primary_roots]

    for lateral_root in lateral_roots
        found_exact_match = false
        for (i, primary_root) in enumerate(primary_roots)

            pr_nb_vertices = reduce(vcat, V.([sg], hypervertex.(V(primary_root))))
            roots_match = [
                sv_lat in pr_nb_vertices && !isnothing(findfirst(se -> src(se) == -3, edges(sv_lat))) &&
                    imag(
                        se_classification_dict[
                            edges(sv_lat)[findfirst(se -> src(se) == -3, edges(sv_lat))],
                        ]
                    ) == 1
                    for sv_lat in V(lateral_root)[[1, end]]
            ] |> all

            if roots_match
                push!(root_systems[i], lateral_root)
                found_exact_match = true
                break
            end
        end

        if !found_exact_match
            i = findmin(pr -> distance(V(lateral_root)[1], pr), primary_roots)[2]
            push!(root_systems[i], lateral_root)
        end
    end

    for rs in root_systems
        sort_root_system!(rs)
    end

    return root_systems
end

function sort_root_system!(rs::Vector{<:Root})
    sort!(rs, by = is_primary, rev = true)
    @assert !is_primary(rs[2]) "Root systems should only contain one primary root"
    sort!(@view(rs[2:end]), by = curve_length, rev = true)

    return nothing
end
