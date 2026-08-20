# # Type
"""
    Root

Represents a single root in a root system.

Root systems, as acquired by [`get_roots`](@ref) are represented by an array of `Root`s.
You can inspect these using the function [`examine`](@ref) or manually calculate the root lengths using [`curve_length`](@ref).
"""
struct Root{T, U}
    is_primary::Bool
    V::Vector{RootVertex{T, U}}
end
function Root(
        re::RootEdge{T, U}, edge_classification_dict::Dict{RootEdge{T, U}, <:Complex},
        rg::RootGraph{T, U}
    ) where {T, U}
    @assert haskey(edge_classification_dict, re)
    return Root(
        real(edge_classification_dict[re]) == 1,
        [getvertex(rg, v) for v in vertices(re)],
    )
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
curve_length(r::Root) = sqrt.(diff(xs(r)) .^ 2 + diff(ys(r)) .^ 2) |> sum
distance(r::Root) = sqrt((ys(r)[end] - ys(r)[1])^2 + (xs(r)[end] - xs(r)[1])^2)
distance(rv::RootVertex, r::Root) = minimum(sqrt.((x(rv) .- xs(r)) .^ 2 + (y(rv) .- ys(r)) .^ 2))
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
    println("Average lateral root length: $(sum(curve_length.(rs[2:end])) / length(rs[2:end]) |> f)")
    println("Individual lateral root lengths:")
    for (i, r) in enumerate(rs)
        i == 1 && continue # only do lateral roots but start from i=2
        println("Root $i: $(curve_length(rs[i]) |> f)")
    end
    return
end

function examine(rss::Vector{<:Vector{<:Root}}; digits = 2)
    for (i, rs) in enumerate(rss)
        println("Root system $i")
        examine(rs; digits)
        println("")
    end
    return
end
