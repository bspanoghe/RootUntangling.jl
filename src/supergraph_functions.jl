polarity(se::SingularEdge{T, U}, sv::SingularVertex{T, U}) where {T, U} = src(se) == id(sv) ? 1 : -1


angle(hv1::HyperVertex{T, U}, hv2::HyperVertex{T, U}; reverse_order::Bool = false) where {T, U} = (
    reverse_order ? atan(y(hv1) - y(hv2), x(hv1) - x(hv2)) :  atan(y(hv2) - y(hv1), x(hv2) - x(hv1))
)
angle(sg::SuperGraph{T, U}, se::SingularEdge{T, U}; reverse_order::Bool = false) where {T, U} = (
    vertices(se) .|> (v -> getsingularvertex(sg, v)) .|> hypervertex |> hvs -> angle(hvs...; reverse_order)
)
angle(sg::SuperGraph{T, U}, he::HyperEdge{T, U}; reverse_order::Bool = false) where {T, U} = (
    vertices(he) .|> (v -> gethypervertex(sg, v)) |> hvs -> angle(hvs...; reverse_order)
)

cosine_similarity(sg::SuperGraph{T, U}, ae1::AbstractEdge{T}, ae2::AbstractEdge{T}, v::T) where {T, U} = (
    cos( -((angle(sg, e, reverse_order = (v == vertices(e)[1])) for e in [ae1, ae2])...) ) # ensure angle of edge is calculated according to same common vertex as starting point
)
cosine_similarity(sg::SuperGraph{T, U}, ae::AbstractEdge{T}, α; reverse_order::Bool = false) where {T, U} = (
    cos( angle(sg, ae; reverse_order) - α )
)

angle_dissimilarity(sg::SuperGraph{T, U}, ae1::AbstractEdge{T}, ae2::AbstractEdge{T}, v::T) where {T, U} = (
    (1 - cosine_similarity(sg, ae1, ae2, v)) / 2
)

xs(sg::SuperGraph{T, U}, se::SingularEdge{T, U}) where {T, U} = x.(V(sg, se))
ys(sg::SuperGraph{T, U}, se::SingularEdge{T, U}) where {T, U} = y.(V(sg, se))