# # edges
"""
    RootEdge{T, U}

Represents the possible presence of a single root along a segment from a root scan.
"""
struct RootEdge{T, U}
    src::T
    dst::T
    width::U
    pred_primary::Union{U, Missing} #! only for testing first 7 ROI - remove Missing Union later

    RootEdge(src::T, dst::T, width::U, pred_primary::Union{U, Missing}) where {T, U} = (
        new{T, U}(minmax(src, dst)..., width, pred_primary)
    )
end
RootEdge(src::T, dst::T) where {T} = RootEdge(src, dst, NaN, NaN)

src(re::RootEdge) = re.src
dst(re::RootEdge) = re.dst
vertices(re::RootEdge) = (src(re), dst(re))
width(re::RootEdge) = re.width
pred_primary(re::RootEdge) = re.pred_primary

is_augmented(v::Integer) = v < 0
is_augmented(re::RootEdge) = any(is_augmented, vertices(re))

# an edge is identified by its vertices
Base.:(==)(re1::RootEdge, re2::RootEdge) = vertices(re1) == vertices(re2)
Base.hash(re::RootEdge, h::UInt) = hash(vertices(re), h)

Base.show(io::IO, re::RootEdge) = print(io, "RootEdge$(vertices(re))")

# # vertices
"""
    RootVertex{T, U}

Represents the possible passage of a single root through a branchpoint or endpoint of a segment in the image.
"""
struct RootVertex{T, U}
    id::T
    edges::Vector{RootEdge{T, U}}
    x::U
    y::U
    pred_split::Union{U, Missing} #! only for testing first 7 ROI - remove later

    RootVertex(id::T, edges::Vector{RootEdge{T, U}}, x::U, y::U, pred_split::Union{U, Missing}) where {T, U} = (
        new{T, U}(id, sort(edges, by = src), x, y, pred_split)
    )
end

id(rv::RootVertex) = rv.id
edges(rv::RootVertex) = rv.edges
edges(rvs::Vector{RootVertex{T, U}}) where {T, U} = reduce(vcat, edges.(rvs), init = RootEdge{T, U}[])
x(rv::RootVertex) = rv.x
y(rv::RootVertex) = rv.y
pred_split(rv::RootVertex) = rv.pred_split

coords(rv::RootVertex) = (x(rv), y(rv))
is_augmented(rv::RootVertex) = is_augmented(id(rv))

Base.show(io::IO, rv::RootVertex) = print(io, "RootVertex$((id(rv), edges(rv)))")

# # segments
# a segment from the scan is represented by the vector of parallel edges lying on it
width(seg::Vector{<:RootEdge}) = width(first(seg))
pred_primary(seg::Vector{<:RootEdge}) = pred_primary(first(seg))
is_augmented(seg::Vector{<:RootEdge}) = is_augmented(first(seg))
srcs(seg::Vector{<:RootEdge}) = unique(src.(seg))
dsts(seg::Vector{<:RootEdge}) = unique(dst.(seg))
vertices(seg::Vector{<:RootEdge}) = [srcs(seg); dsts(seg)]

# # graphs
"""
    RootGraph{T, U}

Represents one or more root systems.

The edges E represent the possible presence of a single root, while the segments S are groups
of edges at the same position, representing a segment in the scan which may contain multiple roots.
Vertices, edges and segments are split into original (V₀ / E₀ / S₀) and augmented (V₊ / E₊ / S₊) ones.
"""
struct RootGraph{T, U}
    V₀::Vector{RootVertex{T, U}}
    V₊::Vector{RootVertex{T, U}}

    E₀::Vector{RootEdge{T, U}}
    E₊::Vector{RootEdge{T, U}}

    S₀::Vector{Vector{RootEdge{T, U}}}
    S₊::Vector{Vector{RootEdge{T, U}}}

    function RootGraph(V₀::Vector{RootVertex{T, U}}, V₊::Vector{RootVertex{T, U}}) where {T, U}
        # sort all vertices
        sort!(V₀, by = id)
        sort!(V₊, by = rv -> -id(rv)) # augmented vertices use negative integers as id

        # assert all ids are consecutive integers starting from 1 / -1
        @assert id.(V₀) == 1:length(V₀)
        @assert id.(V₊) == -1:-1:-length(V₊)

        # get edges and sort
        E₀ = [e for e in edges(V₀) if !is_augmented(e)] |> unique |> es -> sort(es, by = src)
        E₊ = edges(V₊) |> unique |> es -> sort(es, by = src)

        # group edges by position: edges at the same position lie on the same segment of the scan
        position = Dict(id(rv) => coords(rv) for rv in [V₀; V₊])
        S₀ = group(e -> (position[src(e)], position[dst(e)]), E₀)
        S₊ = group(e -> (position[src(e)], position[dst(e)]), E₊)

        return new{T, U}(V₀, V₊, E₀, E₊, S₀, S₊)
    end
end

# group the elements of `xs` by the value of `f`, preserving order
function group(f, xs::Vector{T}) where {T}
    groups = Vector{T}[]
    group_idxs = Dict{Any, Int}()
    for x in xs
        i = get!(() -> (push!(groups, T[]); length(groups)), group_idxs, f(x))
        push!(groups[i], x)
    end
    return groups
end

V₀(rg::RootGraph) = rg.V₀
V₊(rg::RootGraph) = rg.V₊
V(rg::RootGraph) = [V₊(rg); V₀(rg)]

E₀(rg::RootGraph) = rg.E₀
E₊(rg::RootGraph) = rg.E₊
E(rg::RootGraph) = [E₊(rg); E₀(rg)]

segments₀(rg::RootGraph) = rg.S₀
segments₊(rg::RootGraph) = rg.S₊
segments(rg::RootGraph) = [segments₊(rg); segments₀(rg)]

# pairs of edges meeting in a vertex, representing a possible root passing through it
E₂(rv::RootVertex) = [
    [edges(rv)[i], edges(rv)[j]]
        for i in eachindex(edges(rv)) for j in eachindex(edges(rv))
        if (i > j) && !all(is_augmented.(edges(rv)[[i, j]]))
]
E₂(rg::RootGraph) = E₂.(V₀(rg)) |> x -> reduce(vcat, x, init = eltype(x)[])
