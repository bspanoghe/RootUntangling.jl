# get the exclusion sets of a hypervertex `hv`
function exclusion_sets(hv::HyperVertex)
    n_h = num_hypotheses(hv)

    S_idxs = [
        [ones(Int64, level); zeros(Int64, n_h - level)]
        for level in 1:n_h
    ] |> 
        x -> reduce(vcat, x) |>
        cumsum |>
        x -> reshape(x, (n_h, n_h)) |>
        eachrow
    
    S = [vertices(hv)[S_l_idxs] for S_l_idxs in S_idxs]
    
    return S
end

num_hypotheses(hv::HyperVertex) = length(vertices(hv)) |> N -> (sqrt(1 + 8*N) - 1)/2 |> Int64 # from N = (n+1)*n / 2


"""
    order(hv::HyperVertex{T, U}, v::T)

Calculate the order of a vertex `v`. 

Considering the difference in id from every vertex from the hypervertex's root vertex, the lower bound of the set of vertices with order n is sum(1:(n-1)) and therefore the order in function of the difference in id Δx equals floor((1 + sqrt(8* Δx + 1))/2)
"""
order(hv::HyperVertex{T, U}, v::T) where {T, U} = (
    v - vertices(hv)[1] |> x -> (1 + sqrt(8*x + 1))/2 |> x -> floor(Int64, x)
)
order(sv::SingularVertex) = order(hypervertex(sv), id(sv))

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