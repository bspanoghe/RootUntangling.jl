curve_length(r::Root) = sqrt.(diff(xs(r)).^2 + diff(ys(r)).^2) |> sum
distance(r::Root) = sqrt( (ys(r)[end] - ys(r)[1])^2 +  (xs(r)[end] - xs(r)[1])^2 )
tortuosity(r::Root) = curve_length(r) / distance(r)

function switch!(r1::Root, r2::Root, he::HyperEdge)
    # find vertices in roots that match source of hyperedge
    src_idx1 = findfirst(v -> v in vertices(gethypervertex(src(he), Vₕ(sg))), vertices(r1))
    src_idx2 = findfirst(v -> v in vertices(gethypervertex(src(he), Vₕ(sg))), vertices(r2))

    # find vertices in roots that match destination of hyperedge
    dst_idx1 = findfirst(
        v -> v in V(sg, he)[2],
        [get(vertices(r1), i, missing) for i in idx1 .+ [-1, 1]] # only need to look at neighbouring vertices
    )
    dst_idx2 = findfirst(
        v -> v in V(sg, he)[1],
        [get(vertices(r2), i, missing) for i in idx2 .+ [-1, 1]]
    )

    # skip switching if hyperedge is at an extremity of either root (switching does nothing)
    any([idx in [1, length(r1)] for idx in [src_idx1, dst_idx1]]) ||
        any([idx in [1, length(r2)] for idx in [src_idx2, dst_idx2]]) &&
        (@debug("`switch!` skipped"); return nothing)
    
    # check if root orientations match
    orientation_match = (dst_idx1 - src_idx1) == (dst_idx2 - src_idx2)

    # splice roots at hyperedge and switch a half of both
    if orientation_match
        tail1 = splice!(r1.V, dst_idx1:length(r1), r2.V[dst_idx2:end])
        splice!(r2.V, dst_idx2:length(V(r2)), tail1)
    else
        tail1 = splice!(r1.V, dst_idx1:length(r1), reverse(r2.V[1:dst_idx2]))
        splice!(r2.V, 1:dst_idx2, reverse(tail1))        
    end

    return nothing
end

function find_overlaps(sg::SuperGraph, model::JuMP.Model, roots)
    he_classification_dict = get_he_classification_dict(sg, model)
    overlap_hes = filter(he -> imag(he_classification_dict[he]) > 1, Eₕ₀(sg)) #! only looks for lateral overlaps

    overlap_dict = [
        he => [r for r in roots if !any([isdisjoint(vs, vertices(r)) for vs in V(sg, he)])]
        for he in overlap_hes
    ] |> Dict

    return overlap_dict #!

    # if you have multiple edges in a connected linear path, each with the same amount of roots, only choose one edge from that path
    # reasoning: switching roots on multiple of these consecutive edges has no added effect over switching them on just one
    overlap_hes_subset = [
        get_linear_chains(filter(he -> length(overlap_dict[he]) == n, overlap_hes))
        for n in unique(length.(values(overlap_dict)))
    ] |> x -> reduce(vcat, x) .|> first

    overlap_dict_subset = [
        he => [r for r in roots if !any([isdisjoint(vs, vertices(r)) for vs in V(sg, he)])]
        for he in overlap_hes_subset
    ] |> Dict
    
    return overlap_dict_subset
end

# group edges into linear chains
function get_linear_chains(hes::Vector{<:HyperEdge})
    chains = Vector{HyperEdge}[hes[1:1]]
    remaining_hes = hes[2:end]

    while !isempty(remaining_hes)
        chain = chains[end]
        head_idx = findfirst(he -> !isdisjoint(vertices(chain[1]), vertices(he)), remaining_hes)
        if !isnothing(head_idx)
            he = popat!(remaining_hes, head_idx)
            pushfirst!(chain, he)
            continue
        end

        tail_idx = findfirst(he -> !isdisjoint(vertices(chain[end]), vertices(he)), remaining_hes)
        if !isnothing(tail_idx)
            he = popat!(remaining_hes, tail_idx)
            pushfirst!(chain, he)
            
            continue
        end

        he = pop!(remaining_hes)
        push!(chains, [he])
    end

    return chains
end