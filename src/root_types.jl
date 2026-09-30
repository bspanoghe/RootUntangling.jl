# # Types

struct RootFragment{T}
    is_primary::Bool
    edges::Vector{RootArc{T}}
end

is_primary(rf::RootFragment) = rf.is_primary
edges(rf::RootFragment) = rf.edges
is_fullgrown(rf::RootFragment) = all(is_augmented.(edges(rf)[[1, end]]))
vertices(rf::RootFragment) = unique(Iterators.flatten(vertices.(edges(rf))))


"""
    Root

Represents a single root in a root system.
"""
abstract type Root end
is_primary(::Root) = r.is_primary
fragments(::Root) = error("What do you mean you didn't implement the method yet?!")
edges(::Root) = error("What do you mean you didn't implement the method yet?!")
is_fullgrown(::Root) =  error("What do you mean you didn't implement the method yet?!")
vertices(::Root) =  error("What do you mean you didn't implement the method yet?!")
Base.length(r::Root) = length(edges(r))

V(rg::RootGraph, r::Root) = V.([rg], vertices(r))
V₀(rg::RootGraph, r::Root) = V(rg, r)[2:(end-1)] # first and last vertex are always augmented
xs(rg::RootGraph, r::Root) = x.(V₀(rg, r))
ys(rg::RootGraph, r::Root) = y.(V₀(rg, r))
distance(rg::RootGraph, r::Root) = sqrt((ys(rg, r)[end] - ys(rg, r)[1])^2 + (xs(rg, r)[end] - xs(rg, r)[1])^2)
distance(rg::RootGraph, rv::RootVertex, r::Root) = (
    minimum(sqrt.((x(rv) .- xs(rg, r)) .^ 2 + (y(rv) .- ys(rg, r)) .^ 2))
)
curve_length(rg::RootGraph, r::Root) = sqrt.(diff(xs(rg, r)) .^ 2 + diff(ys(rg, r)) .^ 2) |> sum
tortuosity(rg::RootGraph, r::Root) = curve_length(rg, r) / distance(rg, r)
tortuosity(rg::RootGraph, rs::Vector{<:Root}) = sum(tortuosity.([rg], rs))
weighted_tortuosity(rg::RootGraph, r::Root) = curve_length(rg, r) * tortuosity(rg, r)
weighted_tortuosity(rg::RootGraph, rs::Vector{<:Root}) = sum(weighted_tortuosity.([rg], rs))

# root consisting of a single fragment
struct SimpleRoot{T} <: Root
    is_primary::Bool
    fragment::RootFragment{T}
end
fragment(sr::SimpleRoot) = sr.fragment
fragments(sr::SimpleRoot) = [fragment(sr)]
edges(sr::SimpleRoot) = edges(fragment(sr))
is_fullgrown(::SimpleRoot) = true
vertices(sr::SimpleRoot) = vertices(fragment(sr))

struct CompositeRoot{T} <: Root
    is_primary::Bool
    fragments::Vector{RootFragment{T}}
end
fragments(cr::CompositeRoot) = cr.fragments
edges(cr::CompositeRoot) = reduce(vcat, edges.(fragments(cr)))
is_fullgrown(cr::CompositeRoot) = all(is_augmented.( edges(cr)[[1, end]] ))
vertices(cr::CompositeRoot) = unique(reduce(vcat, vertices.(fragments(cr))))

"""
    RootSystem

Represents a complete root system of a single plant. 
See also [`get_rootsystems`](@ref), [`examine`](@ref) and [`curve_length`](@ref).
"""
struct RootSystem{T}
    primary::Union{SimpleRoot{T}, CompositeRoot{T}}
    laterals::Vector{Union{SimpleRoot{T}, CompositeRoot{T}}}
end
primary(rs::RootSystem) = rs.primary
laterals(rs::RootSystem) = rs.laterals
roots(rs::RootSystem) = [primary(rs); laterals(rs)]
Base.length(rs::RootSystem) = length(roots(rs))

tortuosity(rg::RootGraph, rs::RootSystem) = tortuosity(rg, primary(rs)) + tortuosity(rg, laterals(rs))
tortuosity(rg::RootGraph, rss::Vector{<:RootSystem}) = sum(tortuosity.([rg], rss))
weighted_tortuosity(rg::RootGraph, rs::RootSystem) = weighted_tortuosity(rg, primary(rs)) + weighted_tortuosity(rg, laterals(rs))
weighted_tortuosity(rg::RootGraph, rss::Vector{<:RootSystem}) = sum(weighted_tortuosity.([rg], rss))

"""
    examine(rg, rs)

Get functional information from root systems.
"""
function examine end

function examine(rg::RootGraph, rs::RootSystem; digits = 2)
    f = x -> round(x; digits)

    println("Primary root length: $(curve_length(rg, primary(rs)) |> f)")
    println("Number of lateral roots: $(length(laterals(rs)))")
    println("Average lateral root length: $(sum(curve_length.([rg], laterals(rs))) / length(laterals(rs)) |> f)")
    println("Individual lateral root lengths:")
    for (i, r) in enumerate(laterals(rs))
        println("Root $i: $(curve_length(rg, r) |> f)")
    end
    return
end

function examine(rg::RootGraph, rss::Vector{<:RootSystem}; digits = 2)
    for (i, rs) in enumerate(rss)
        println("Root system $i")
        examine(rg, rs; digits)
        println("")
    end
    return
end