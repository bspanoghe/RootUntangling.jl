# # Type
struct RootFragment{T, U}
    is_primary::Bool
    edges::Vector{RootEdge{T, U}}
    start::RootVertex{T, U}
    stop::RootVertex{T, U}
end
is_primary(rf::RootFragment) = rf.is_primary
edges(rf::RootFragment) = rf.edges
start(rf::RootFragment) = rf.start
stop(rf::RootFragment) = rf.stop

"""
    Root

Represents a single root in a root system.
"""
struct Root{T, U}
    is_primary::Bool
    edges::Vector{RootEdge{T, U}}
    correct_directions::Vector{Bool}
    # V::Vector{RootVertex{T, U}}
end
is_primary(r::Root) = r.is_primary
edges(r::Root) = r.edges
correct_directions(r::Root) = r.correct_directions

function V(r::Root)
    es = edges(r)
    correct_dir = correct_directions(r)
    return [
        correct_dir[i] ? V(es[i]) : reverse(V(es[i])) 
        for i in eachindex(es)
    ] |> unique
end

vertices(r::Root) = id.(V(r))
xs(r::Root) = x.(V(r))
ys(r::Root) = y.(V(r))
function Eₕ(r::Root)
    root_rvs = rootvertex.(V(r))
    return [RootEdge(id.(root_rvs)[i], id.(root_rvs)[i-1]) for i in 2:length(root_rvs)]
end

"""
    RootSystem

Represents a complete root system of a single plant. 
See also [`get_rootsystems`](@ref), [`examine`](@ref) and [`curve_length`](@ref).
"""
struct RootSystem{T, U}
    primary::Root{T, U}
    laterals::Vector{Root{T, U}}
end
primary(rs::RootSystem) = rs.primary
laterals(rs::RootSystem) = rs.laterals
roots(rs::RootSystem) = [primary(rs); laterals(rs)]

"""
    curve_length(r::Root)

Calculate the physical length of a root.
"""
curve_length(r::Root) = sqrt.(diff(xs(r)) .^ 2 + diff(ys(r)) .^ 2) |> sum
distance(r::Root) = sqrt((ys(r)[end] - ys(r)[1])^2 + (xs(r)[end] - xs(r)[1])^2)
distance(v::RootVertex, r::Root) = minimum(sqrt.((x(v) .- xs(r)) .^ 2 + (y(v) .- ys(r)) .^ 2))
tortuosity(r::Root) = curve_length(r) / distance(r)
tortuosity(rs::RootSystem) = sum(tortuosity.(roots(rs)))

function roughness(r::Root)
    angles = [angle(V(r)[i], V(r)[i-1]) for i in 2:length(r)]
    length(angles) == 1 && (return 0.0)
    angle_changes = diff(angles)
    return sum(angle_changes.^2) / length(angle_changes)
end
roughness(rs::RootSystem) = sum(roughness.(roots(rs)))


"""
    examine(rs)

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