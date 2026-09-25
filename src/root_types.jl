# # Types
abstract type RootFragment end
is_primary(rf::RootFragment) = rf.is_primary
edges(rf::RootFragment) = rf.edges
is_fullgrown(rf::RootFragment) = all(is_augmented.(edges(rf)[[1, end]]))
vertices(rf::RootFragment) = unique(Iterators.flatten(vertices.(edges(rf))))

mutable struct DirectedRootFragment{T} <: RootFragment
    is_primary::Bool
    edges::Vector{RootArc{T}}
end

mutable struct UndirectedRootFragment{T, U} <: RootFragment
    is_primary::Bool
    edges::Vector{RootEdge{T, U}}
end

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

# root consisting of a single fragment
struct SimpleRoot{T} <: Root
    is_primary::Bool
    fragment::DirectedRootFragment{T}
end
fragment(sr::SimpleRoot) = sr.fragment
fragments(sr::SimpleRoot) = [fragment(sr)]
edges(sr::SimpleRoot) = edges(fragment(sr))
is_fullgrown(::SimpleRoot) = true
vertices(sr::SimpleRoot) = vertices(fragment(sr))

struct CompositeRoot <: Root
    is_primary::Bool
    fragments::Vector{<:RootFragment}
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
struct RootSystem
    primary::Root
    laterals::Vector{<:Root}
    function RootSystem(primary::Root, laterals::Vector{<:Root})
        return new(primary, sort(laterals, by = length, rev = true))
    end
end
primary(rs::RootSystem) = rs.primary
laterals(rs::RootSystem) = rs.laterals
roots(rs::RootSystem) = [primary(rs); laterals(rs)]
Base.length(rs::RootSystem) = length(roots(rs))


# tortuosity(rs::RootSystem) = sum(tortuosity.(roots(rs)))

# function roughness(r::Root)
#     angles = [angle(V(r)[i], V(r)[i-1]) for i in 2:length(r)]
#     length(angles) == 1 && (return 0.0)
#     angle_changes = diff(angles)
#     return sum(angle_changes.^2) / length(angle_changes)
# end
# roughness(rs::RootSystem) = sum(roughness.(roots(rs)))


# """
#     examine(rs)

# Get functional information from root systems.
# """
# function examine end

# function examine(rs::RootSystem; digits = 2)
#     f = x -> round(x; digits)

#     println("Primary root length: $(curve_length(primary(rs)) |> f)")
#     println("Number of lateral roots: $(length(laterals(rs)))")
#     println("Average lateral root length: $(sum(curve_length.(laterals(rs))) / length(laterals(rs)) |> f)")
#     println("Individual lateral root lengths:")
#     for (i, r) in enumerate(laterals(rs))
#         println("Root $i: $(curve_length(r) |> f)")
#     end
#     return
# end

# function examine(rss::Vector{<:RootSystem}; digits = 2)
#     for (i, rs) in enumerate(rss)
#         println("Root system $i")
#         examine(rs; digits)
#         println("")
#     end
#     return
# end