Base.length(r::Root) = sqrt.(diff(xs(r)).^2 + diff(ys(r)).^2) |> sum
distance(r::Root) = sqrt( (ys(r)[end] - ys(r)[1])^2 +  (xs(r)[end] - xs(r)[1])^2 )
tortuosity(r::Root) = length(r) / distance(r)

function switch!(r1::Root, r2::Root, v1::T, v2::T) where {T}
    idx1 = findfirst(v -> v == v1, vertex_ids(r1))
    idx2 = findfirst(v -> v == v1, vertex_ids(r1))
end