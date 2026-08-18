# get singular / hypervertex from vertex
getsingularvertex(sg::SuperGraph{T, U}, v::T) where {T, U} = v > 0 ? V₀(sg)[v] : V₊(sg)[-v]
gethypervertex(Vₕ::Vector{HyperVertex{T, U}}, v::T) where {T, U} = Vₕ[findfirst(hv -> id(hv) == v, Vₕ)]
gethypervertex(sg::SuperGraph{T, U}, v::T) where {T, U} = v > 0 ? Vₕ₀(sg)[v] : Vₕ₊(sg)[-v]

# get information from a singularedge that requires graph information
xs(sg::SuperGraph{T, U}, se::SingularEdge{T, U}) where {T, U} = x.(V(sg, se))
ys(sg::SuperGraph{T, U}, se::SingularEdge{T, U}) where {T, U} = y.(V(sg, se))
width(sg::SuperGraph, se::SingularEdge) = width(hyperedge(se)) / minimum(order.(getsingularvertex.([sg], vertices(se))))

# get neighbors
neighbor(hv::HyperVertex{T, U}, he::HyperEdge{T, U}, Vₕ::Vector{HyperVertex{T, U}}) where {T, U} = (
    vertices(he)[findfirst(v -> v != id(hv), vertices(he))] |> (v -> gethypervertex(Vₕ, v))
)
neighbor(sg::SuperGraph{T, U}, sv::SingularVertex{T, U}, se::SingularEdge{T}) where {T, U} = (
    vertices(se)[findfirst(v -> v != id(sv), vertices(se))] |> v -> getsingularvertex(sg, v)
)
neighbor(sg::SuperGraph{T, U}, hv::HyperVertex{T, U}, he::HyperEdge{T}) where {T, U} = (
    vertices(he)[findfirst(v -> v != id(hv), vertices(he))] |> v -> gethypervertex(sg, v)
)
neighbors(sg::SuperGraph{T, U}, av::AbstractVertex{T, U}) where {T, U} = [neighbor(sg, av, ae) for ae in edges(av)]

inner_vertices(sg::SuperGraph) = [v for v in V₀(sg) if length([n for n in neighbors(sg, Vₕ(v)) if !is_augmented(n)]) > 1]
outer_vertices(sg::SuperGraph) = [v for v in V₀(sg) if length([n for n in neighbors(sg, Vₕ(v)) if !is_augmented(n)]) == 1]

# polarity
polarity(se::SingularEdge{T, U}, sv::SingularVertex{T, U}) where {T, U} = src(se) == id(sv) ? 1 : -1
polarity(sv1::SingularVertex{T, U}, sv2::SingularVertex{T, U}) where {T, U} = id(sv1) < id(sv2) ? 1 : -1

# angles
angle(hv1::HyperVertex{T, U}, hv2::HyperVertex{T, U}; reverse_order::Bool = false) where {T, U} = (
    reverse_order ? atan(y(hv1) - y(hv2), x(hv1) - x(hv2)) : atan(y(hv2) - y(hv1), x(hv2) - x(hv1))
)
angle(sg::SuperGraph{T, U}, se::SingularEdge{T, U}; reverse_order::Bool = false) where {T, U} = (
    vertices(se) .|> (v -> getsingularvertex(sg, v)) .|> hypervertex |> hvs -> angle(hvs...; reverse_order)
)
angle(sg::SuperGraph{T, U}, he::HyperEdge{T, U}; reverse_order::Bool = false) where {T, U} = (
    vertices(he) .|> (v -> gethypervertex(sg, v)) |> hvs -> angle(hvs...; reverse_order)
)

cosine_similarity(sg::SuperGraph{T, U}, ae1::AbstractEdge{T}, ae2::AbstractEdge{T}, v::T) where {T, U} = (
    cos(-((angle(sg, e, reverse_order = (v == vertices(e)[1])) for e in [ae1, ae2])...)) # ensure angle of edge is calculated according to same common vertex as starting point
)
cosine_similarity(sg::SuperGraph{T, U}, ae::AbstractEdge{T}, α; reverse_order::Bool = false) where {T, U} = (
    cos(angle(sg, ae; reverse_order) - α)
)

angle_dissimilarity(sg::SuperGraph{T, U}, ae1::AbstractEdge{T}, ae2::AbstractEdge{T}, v::T) where {T, U} = (
    (1 - cosine_similarity(sg, ae1, ae2, v)) / 2
)
