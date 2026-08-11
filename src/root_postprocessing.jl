Base.length(r::Root) = sqrt.(diff(xs(r)).^2 + diff(ys(r)).^2) |> sum
distance(r::Root) = sqrt( (ys(r)[end] - ys(r)[1])^2 +  (xs(r)[end] - xs(r)[1])^2 )
tortuosity(r::Root) = length(r) / distance(r)

function switch!(sg::SuperGraph, r1::Root, r2::Root, v1::T, v2::T) where {T}
    idx1 = findfirst(v -> v == v1, vertex_ids(r1))
    idx2 = findfirst(v -> v == v2, vertex_ids(r2))

    v1_downstream_nb = vertex_ids(r1)[idx1+1]
    v2_downstream_nb = vertex_ids(r2)[idx2+1]

    orientation_match = isequal(
        hypervertex(getsingularvertex(sg, v1_downstream_nb)),
        hypervertex(getsingularvertex(sg, v2_downstream_nb))
    )

    if orientation_match
        
    else

    end

    return nothing
end