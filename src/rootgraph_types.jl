# # ids
is_augmented(v::Integer) = v < 0

# # edges
"""
    RootEdge{T, U}

Represents a segment from a root scan, which may contain roots going in either direction.
"""
struct RootEdge{T, U}
    src::T
    dst::T
    segment_id::T
    width::U
    pred_primary::Union{U, Missing} #! only for testing first 7 ROI - remove Missing Union later

    RootEdge(src::T, dst::T, segment_id::T, width::U, pred_primary::Union{U, Missing}) where {T, U} = (
        new{T, U}(sort([src, dst])..., segment_id, width, pred_primary)
    )
end
RootEdge(src::T, dst::T) where {T} = RootEdge(src, dst, 0, NaN, NaN)
RootEdge(s::Segment) = RootEdge(vertices(s)..., id(s), width(s), pred_primary(s))

src(re::RootEdge) = re.src
dst(re::RootEdge) = re.dst
segment_id(re::RootEdge) = re.segment_id
width(re::RootEdge) = re.width
pred_primary(re::RootEdge) = re.pred_primary
vertices(re::RootEdge) = (src(re), dst(re))

is_augmented(re::RootEdge) = any(is_augmented.(vertices(re)))

Base.sort(res::Vector{<:RootEdge}) = sort(res, by = vertices)
Base.:(==)(re1::RootEdge, re2::RootEdge) = vertices(re1) == vertices(re2)
Base.unique(res::Vector{<:RootEdge}) = unique(x -> vertices(x), res) # doesn't automatically use my equality operator :(
Base.unique!(res::Vector{<:RootEdge}) = unique!(x -> vertices(x), res) # sad times
Base.show(io::IO, re::RootEdge) = print(io, "$(typeof(re).name.name)$(vertices(re))")
Base.show(io::IO, res::Vector{<:RootEdge}) = print(io, "$(typeof(res).name.name)$(vertices.(res))")

# vertices
"""
    RootVertex{T, U}

Represents a branchpoint or endpoint of a segment in the image.
"""
struct RootVertex{T, U}
    id::T
    edges::Vector{RootEdge{T, U}}
    x::U
    y::U
    pred_split::Union{U, Missing} #! only for testing first 7 ROI - remove later

    RootVertex(id::T, edges::Vector{RootEdge{T, U}}, x::U, y::U, pred_split::Union{U, Missing}) where {T, U} = (
        new{T, U}(id, sort(edges, by = e -> src(e)), x, y, pred_split)
    )
end
RootVertex(id, edges, x, y, pred_split) = RootVertex(id, edges, x, y, pred_split)

id(rv::RootVertex) = rv.id
edges(rv::RootVertex) = rv.edges
x(rv::RootVertex) = rv.x
y(rv::RootVertex) = rv.y
pred_split(rv::RootVertex) = rv.pred_split

edges(rvs::Vector{<:RootVertex}) = reduce(vcat, edges.(rvs), init = eltype(rvs)[])
coords(rv::RootVertex) = (x(rv), y(rv))
is_augmented(rv::RootVertex) = is_augmented(id(rv))

Base.show(io::IO, rv::RootVertex) = print(io, "$(typeof(rv).name.name)$((id(rv), edges(rv)))")
Base.show(io::IO, rvs::Vector{<:RootVertex}) = print(io, "$(typeof(rvs).name.name)$(id.(rvs))")

# graphs
"""
    RootGraph{T, U}

Represents one or more root systems.
"""
struct RootGraph{T, U}
    V₀::Vector{RootVertex{T, U}}
    V₊::Vector{RootVertex{T, U}}
    E₀::Vector{RootEdge{T, U}}
    E₊::Vector{RootEdge{T, U}}

    function RootGraph(
            V₀::Vector{RootVertex{T, U}}, V₊::Vector{RootVertex{T, U}}
        ) where {T, U}

        # sort all vertices
        sort!(V₀, by = x -> id(x))
        sort!(V₊, by = x -> -id(x)) # special vertices use negative integers as id

        # assert all ids are consecutive integers starting from 1 / -1
        @assert id.(V₀) == 1:length(V₀)
        @assert id.(V₊) == -1:-1:-length(V₊)

        # get edges and sort
        E₀ = [e for e in edges(V₀) if !(is_augmented(e))] |> unique |> es -> sort(es, by = e -> src(e))
        E₊ = edges(V₊) |> unique |> es -> sort(es, by = e -> src(e))

        return new{T, U}(V₀, V₊, E₀, E₊)
    end
end
V₀(rg::RootGraph) = rg.V₀
V₊(rg::RootGraph) = rg.V₊
E₀(rg::RootGraph) = rg.E₀
E₊(rg::RootGraph) = rg.E₊

V(rg::RootGraph) = [V₊(rg); V₀(rg)]
E(rg::RootGraph) = [E₊(rg); E₀(rg)]

Base.length(rg::RootGraph) = length(V₀(rg))

E₂(rv::RootVertex) = edges(rv) |> es -> [
    sort([es[i], es[j]])
    for i in 1:length(es) for j in i+1:length(es)
    if !all(is_augmented.(es[[i, j]]))
]
E₂(rg::RootGraph) = E₂.(V₀(rg)) |> x -> reduce(vcat, x, init = eltype(x)[])
E₂(rv::RootVertex, re::RootEdge) = edges(rv) |> es -> [
    sort([re, nb_re])
    for nb_re in es
    if re != nb_re && !all(is_augmented.([re, nb_re]))
]
E₂(rg::RootGraph, re::RootEdge) = [
    sort([re, nb_re])
    for nb_re in unique(reduce(vcat, edges.(V(rg, re))))
    if nb_re != re && !all(is_augmented.([re, nb_re]))
]


# additional methods using mathematical syntax of V / V₀ and E / E₀
V(rg::RootGraph{T, U}, re::RootEdge{T, U}) where {T, U} = [getrootvertex(rg, v) for v in vertices(re)]
V(rg::RootGraph{T, U}, rv::RootVertex{T, U}) where {T, U} = [getrootvertex(rg, v) for v in vertices(rv)]

E(rv::RootVertex) = edges(rv)
E(rvs::Vector{<:RootVertex}) = reduce(vcat, E.(rvs), init = eltype(rvs)[])
E(rg::RootGraph{T, U}, rv::RootVertex{T, U}) where {T, U} = E(V(rg, rv))
E₀(rv::RootVertex) = filter(!is_augmented, E(rv))