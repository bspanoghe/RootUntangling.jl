# get a rootvertex from vertex / edge
getvertex(rg::RootGraph{T, U}, v::T) where {T, U} = v > 0 ? V₀(rg)[v] : V₊(rg)[-v]
V(rg::RootGraph{T, U}, v::T) where {T, U} = getvertex(rg, v)
V(rg::RootGraph{T, U}, re::RootEdge{T, U}) where {T, U} = [getvertex(rg, v) for v in vertices(re)]

# get edge coordinates
xs(rg::RootGraph{T, U}, re::RootEdge{T, U}) where {T, U} = x.(V(rg, re))
ys(rg::RootGraph{T, U}, re::RootEdge{T, U}) where {T, U} = y.(V(rg, re))

# get neighbors
neighbor(rg::RootGraph{T, U}, rv::RootVertex{T, U}, re::RootEdge{T, U}) where {T, U} = (
    vertices(re)[findfirst(v -> v != id(rv), vertices(re))] |> v -> getvertex(rg, v)
)
neighbors(rg::RootGraph{T, U}, rv::RootVertex{T, U}) where {T, U} = [neighbor(rg, rv, re) for re in edges(rv)]

inner_vertices(rg::RootGraph) = [rv for rv in V₀(rg) if length(filter(!is_augmented, neighbors(rg, rv))) > 1]
outer_vertices(rg::RootGraph) = [rv for rv in V₀(rg) if length(filter(!is_augmented, neighbors(rg, rv))) == 1]

# polarity
polarity(re::RootEdge{T, U}, rv::RootVertex{T, U}) where {T, U} = src(re) == id(rv) ? 1 : -1
polarity(rv1::RootVertex{T, U}, rv2::RootVertex{T, U}) where {T, U} = id(rv1) < id(rv2) ? 1 : -1

# angles
angle(rv1::RootVertex{T, U}, rv2::RootVertex{T, U}; reverse_order::Bool = false) where {T, U} = (
    reverse_order ? atan(y(rv1) - y(rv2), x(rv1) - x(rv2)) : atan(y(rv2) - y(rv1), x(rv2) - x(rv1))
)
angle(rg::RootGraph{T, U}, re::RootEdge{T, U}; reverse_order::Bool = false) where {T, U} = (
    angle(V(rg, re)...; reverse_order)
)

cosine_similarity(rg::RootGraph{T, U}, re1::RootEdge{T, U}, re2::RootEdge{T, U}, v::T) where {T, U} = (
    cos(-((angle(rg, e, reverse_order = (v == vertices(e)[1])) for e in [re1, re2])...)) # ensure angle of edge is calculated according to same common vertex as starting point
)
cosine_similarity(rg::RootGraph{T, U}, re::RootEdge{T, U}, α; reverse_order::Bool = false) where {T, U} = (
    cos(angle(rg, re; reverse_order) - α)
)

angle_dissimilarity(rg::RootGraph{T, U}, re1::RootEdge{T, U}, re2::RootEdge{T, U}, v::T) where {T, U} = (
    (1 - cosine_similarity(rg, re1, re2, v)) / 2
)
