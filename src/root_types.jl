# # Types

# Like a RootEdge but then directed!
struct RootArc{T, U} <: AbstractEdge
    src::T
    dst::T
    segment_id::T
    xs::Vector{U}
    ys::Vector{U}
end
function RootArc(rg::RootGraph, re::RootEdge, s::Segment; keep_order::Bool)
    source, dest = keep_order ? (src(re), dst(re)) : (dst(re), src(re))
    is_augmented(re) && return RootArc(source, dest, segment_id(re), xs(s), ys(s))
    
    dist_source = distance(x(V(rg, source)), y(V(rg, source)), xs(s)[1], ys(s)[1])
    dist_dest = distance(x(V(rg, dest)), y(V(rg, dest)), xs(s)[1], ys(s)[1])
    xs_cor, ys_cor = (dist_source < dist_dest) ? (xs(s), ys(s)) : (reverse(xs(s)), reverse(ys(s)))

    RootArc(source, dest, segment_id(re), xs_cor, ys_cor)
end
xs(ra::RootArc) = ra.xs
ys(ra::RootArc) = ra.ys

struct RootFragment{T, U}
    is_primary::Bool
    edges::Vector{RootArc{T, U}}
end

is_primary(rf::RootFragment) = rf.is_primary
edges(rf::RootFragment) = rf.edges
is_fullgrown(rf::RootFragment) = all(is_augmented.(edges(rf)[[1, end]]))
vertices(rf::RootFragment) = vertices.(edges(rf)) |> vs -> [first.(vs); last(vs[end])]
xs(rf::RootFragment) = reduce(vcat, xs.(edges(rf)))
ys(rf::RootFragment) = reduce(vcat, ys.(edges(rf)))

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
xs(::Root) =  error("What do you mean you didn't implement the method yet?!")
ys(::Root) =  error("What do you mean you didn't implement the method yet?!")

Base.length(r::Root) = length(edges(r))

distance(x1::T, y1::T, x2::T, y2::T) where {T<:Number} = sqrt((y2-y1)^2 + (x2-x1)^2)
distance(r::Root) = distance(xs(r)[1], ys(r)[1], xs(r)[end], ys(r)[end])
distance(r::Root, x, y) = minimum(distance.([x], xs(r), [y], ys(r)))
curve_length(r::Root) = sqrt.(diff(xs(r)) .^ 2 + diff(ys(r)) .^ 2) |> sum
tortuosity(r::Root) = curve_length(r) / distance(r)
tortuosity(rs::Vector{<:Root}) = sum(tortuosity.(rs))
weighted_tortuosity(r::Root) = curve_length(r) * tortuosity(r)
weighted_tortuosity(rs::Vector{<:Root}) = sum(weighted_tortuosity.(rs))

# root consisting of a single fragment
struct SimpleRoot{T, U} <: Root
    is_primary::Bool
    fragment::RootFragment{T, U}
end
fragment(sr::SimpleRoot) = sr.fragment
fragments(sr::SimpleRoot) = [fragment(sr)]
edges(sr::SimpleRoot) = edges(fragment(sr))
is_fullgrown(::SimpleRoot) = true
vertices(sr::SimpleRoot) = vertices(fragment(sr))
xs(sr::SimpleRoot) = xs(fragment(sr))
ys(sr::SimpleRoot) = ys(fragment(sr))

struct CompositeRoot{T, U} <: Root
    is_primary::Bool
    fragments::Vector{RootFragment{T, U}}
end
fragments(cr::CompositeRoot) = cr.fragments
edges(cr::CompositeRoot) = reduce(vcat, edges.(fragments(cr)))
is_fullgrown(cr::CompositeRoot) = all(is_augmented.( edges(cr)[[1, end]] ))
vertices(cr::CompositeRoot) = unique(reduce(vcat, vertices.(fragments(cr))))
xs(cr::CompositeRoot) = reduce(vcat, xs.([
    edges(fragments(cr)[1]);
    [edges(fragments(cr)[i])[2:end] for i in 2:length(fragments(cr))]... # first edge of subsequent fragment equals last edge of the previous
]))
ys(cr::CompositeRoot) = reduce(vcat, ys.([
    edges(fragments(cr)[1]);
    [edges(fragments(cr)[i])[2:end] for i in 2:length(fragments(cr))]... # first edge of subsequent fragment equals last edge of the previous
]))

"""
    RootSystem

Represents a complete root system of a single plant. 
See also [`get_rootsystems`](@ref), [`examine`](@ref) and [`curve_length`](@ref).
"""
struct RootSystem{T, U}
    primary::Union{SimpleRoot{T, U}, CompositeRoot{T, U}}
    laterals::Vector{Union{SimpleRoot{T, U}, CompositeRoot{T, U}}}
end
primary(rs::RootSystem) = rs.primary
laterals(rs::RootSystem) = rs.laterals
roots(rs::RootSystem) = [primary(rs); laterals(rs)]
Base.length(rs::RootSystem) = length(roots(rs))

tortuosity(rs::RootSystem) = tortuosity(primary(rs)) + tortuosity(laterals(rs))
tortuosity(rss::Vector{<:RootSystem}) = sum(tortuosity.(rss))
weighted_tortuosity(rs::RootSystem) = weighted_tortuosity(primary(rs)) + weighted_tortuosity(laterals(rs))
weighted_tortuosity(rss::Vector{<:RootSystem}) = sum(weighted_tortuosity.(rss))

"""
    examine(rg, rs)

Get functional information from root systems.
"""
function examine end

function examine(rs::RootSystem; digits = 2)
    f = x -> round(x; digits)

    println("Primary root length: $(curve_length(primary(rs)) |> f)")
    println("Number of lateral roots: $(length(laterals(rs)))")
    println("Average lateral root length: $(sum(curve_length.(laterals(rs))) / length(laterals(rs)) |> f)")
    println("Individual lateral root lengths:")
    for (i, r) in enumerate(laterals(rs))
        println("Root $i: $(curve_length(r) |> f)")
    end
    return
end

function examine(rss::Vector{<:RootSystem}; digits = 2)
    for (i, rs) in enumerate(rss)
        println("Root system $i")
        examine(rs; digits)
        println("")
    end
    return
end