# get rootvertex from its id ("the vertex")
getrootvertex(rg::RootGraph{T, U}, v::T) where {T, U} = v > 0 ? V₀(rg)[v] : V₊(rg)[-v]

# get coords from edges
xs(rg::RootGraph{T, U}, re::RootEdge{T, U}) where {T, U} = x.(V(rg, re))
ys(rg::RootGraph{T, U}, re::RootEdge{T, U}) where {T, U} = y.(V(rg, re))

# get neighbors
neighbor(rv::RootVertex{T, U}, re::RootEdge{T, U}) where {T, U} = (
    vertices(re)[findfirst(v -> v != id(rv), vertices(re))]
)
neighbor(rg::RootGraph{T, U}, rv::RootVertex{T, U}, re::RootEdge{T, U}) where {T, U} = (
    getrootvertex(rg, neighbor(rv, re))
)
neighbors(rg::RootGraph{T, U}, rv::RootVertex{T, U}) where {T, U} = [neighbor(rg, rv, re) for re in edges(rv)]

inner_vertices(rg::RootGraph) = [v for v in V₀(rg) if length([n for n in neighbors(rg, v) if !is_augmented(n)]) > 1]
outer_vertices(rg::RootGraph) = [v for v in V₀(rg) if length([n for n in neighbors(rg, v) if !is_augmented(n)]) == 1]

# direction
direction(rv::RootVertex{T, U}, re::RootEdge{T}) where {T, U} = src(re) == id(rv) ? 1 : -1
direction(rv1::RootVertex{T, U}, rv2::RootVertex{T, U}) where {T, U} = id(rv1) < id(rv2) ? 1 : -1
direction(c::Vector{RootEdge{T, U}}, e::RootEdge{T, U}) where {T, U} = e == c[2] ? 1 : -1

# angles
angle(rv1::RootVertex{T, U}, rv2::RootVertex{T, U}) where {T, U} = angle(x(rv1), y(rv1), x(rv2), y(rv2))
angle(rv1::RootVertex{T, U}, rv2::RootVertex{T, U}, correct_order::Bool) where {T, U} = (
    correct_order ? angle(rv1, rv2) : angle(rv2, rv1)
)
angle(rg::RootGraph{T, U}, re::RootEdge{T, U}, correct_order::Bool = true) where {T, U} = (
    angle(V(rg, re)..., correct_order)
)
cosine_similarity(rg::RootGraph{T, U}, re::RootEdge{T, U}, α::Real, correct_order::Bool) where {T, U} = (
    cos(angle(rg, re, correct_order) - α)
)

local_angle(re::RootEdge{T, U}, correct_order::Bool) where {T, U} = (
    correct_order ? src_angle(re) : dst_angle(re)
)
local_cosine_similarity(re1::RootEdge{T, U}, re2::RootEdge{T, U}, v::T) where {T, U} = (
    cos(-((local_angle(re, (v == src(re))) for re in [re1, re2])...)) # ensure angle of edge is calculated according to same common vertex as starting point
)
angle_dissimilarity(re1::RootEdge{T, U}, re2::RootEdge{T, U}, v::T) where {T, U} = (
    (1 - local_cosine_similarity(re1, re2, v)) / 2
)