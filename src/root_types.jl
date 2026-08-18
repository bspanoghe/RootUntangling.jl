# # Type
"""
    Root

Represents a single root in a root system.

Root systems, as acquired by [`get_roots`](@ref) are represented by an array of `Root`s.
You can inspect these using the function [`examine`](@ref) or manually calculate the root lengths using [`curve_length`](@ref).
"""
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
"""
    curve_length(r::Root)

Calculate the physical length of a root.
"""
curve_length(r::Root) = sqrt.(diff(xs(r)).^2 + diff(ys(r)).^2) |> sum
distance(r::Root) = sqrt( (ys(r)[end] - ys(r)[1])^2 +  (xs(r)[end] - xs(r)[1])^2 )
distance(v::SingularVertex, r::Root) = minimum(sqrt.( (x(v) .- xs(r)).^2 + (y(v) .- ys(r)).^2 ))
tortuosity(r::Root) = curve_length(r) / distance(r)
tortuosity(rs::Vector{<:Root}) = sum(tortuosity.(rs))

"""
    examine(rs)

Get functional information from root systems.
"""
function examine end

function examine(rs::Vector{<:Root}; digits = 2)
    @assert length(filter(is_primary, rs)) == 1 "A vector of `Root`s should only have one primary root"
    @assert issorted(rs, by = is_primary, rev = true)

    f = x -> round(x; digits)

    println("Primary root length: $(curve_length(rs[1]) |> f)")
    println("Number of lateral roots: $(length(rs) - 1)")
    println("Average lateral root length: $( sum(curve_length.(rs[2:end]))/length(rs[2:end]) |> f)")
    println("Individual lateral root lengths:")
    for (i, r) in enumerate(rs)
        i == 1 && continue # only do lateral roots but start from i=2
        println("Root $i: $(curve_length(rs[i]) |> f)")
    end
end

function examine(rss::Vector{<:Vector{<:Root}}; digits = 2)
    for (i, rs) in enumerate(rss)
        println("Root system $i")
        examine(rs; digits)
        println("")
    end
end

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
        sort!(roots, by = r -> is_primary(r), rev = true)
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
                sv_lat in pr_nb_vertices && 
                imag(se_classification_dict[
                    edges(sv_lat)[findfirst(se -> src(se) == -3, edges(sv_lat))]
                    ]) == 1
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
        sort!(rs, by = r -> is_primary(r), rev = true)
    end

    return root_systems
end