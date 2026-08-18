function greedy_switch(sg, model, roots; max_tries = 100)
    roots_copy = deepcopy(roots)
    current_tort = tortuosity(roots)
    improving = true

    counter = 0
    while improving
        counter += 1
        counter > max_tries && (@info "Max tries reached"; break)

        overlap_dict = find_overlaps(sg, model, roots_copy)
        switch_dict = get_switch_dict(overlap_dict, roots_copy)
        n = maximum(keys(switch_dict))

        torts = [total_tortuosity(sg, roots_copy, create_switches(n, i), switch_dict) for i in 1:n]
        best_idx = argmin(torts)
        if torts[best_idx] < current_tort
            make_switches!(sg, roots_copy, create_switches(n, best_idx), switch_dict)
        else
            improving = false
        end
    end

    return roots_copy
end

# find hyperedges where multiple roots overlap (and can switch)
function find_overlaps(sg::SuperGraph, model::JuMP.Model, roots)
    he_classification_dict = get_he_classification_dict(sg, model)
    overlap_hes = filter(he -> imag(he_classification_dict[he]) > 1, Eₕ₀(sg)) #! only looks for lateral overlaps

    # map all hyperedges to the roots they are part of
    # discarding roots of length 2 or smaller (switching does nothing)
    overlap_dict = [
        he => [r for r in roots if !any([isdisjoint(vs, vertices(r)) for vs in V(sg, he)]) && length(r) > 2]
            for he in overlap_hes
    ] |> Dict

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

function get_switch_dict(overlap_dict, roots)
    switch_pairs = Pair[]
    counter = 0
    for (he, overlap_roots) in overlap_dict
        for i in eachindex(overlap_roots), j in eachindex(overlap_roots)
            j <= i && continue # order of roots is not important
            counter += 1
            r1_idx = findfirst(r -> r == overlap_roots[i], roots)
            r2_idx = findfirst(r -> r == overlap_roots[j], roots)
            push!(switch_pairs, counter => (he, r1_idx, r2_idx))
        end
    end

    return Dict(switch_pairs)
end

function total_tortuosity(sg, roots, switches, switch_dict)
    roots_copy = deepcopy(roots)
    make_switches!(sg, roots_copy, switches, switch_dict)
    return tortuosity(roots_copy)
end

create_switches(n::Integer, i::Integer) = [zeros(Bool, i - 1); true; zeros(Bool, n - i)]

function make_switches!(sg, roots, switches, switch_dict)
    for switch in findall(switches)
        he, r1_idx, r2_idx = switch_dict[switch]
        switch!(sg, he, roots[r1_idx], roots[r2_idx])
    end

    return nothing
end

# perform a crossing over between two roots
function switch!(sg::SuperGraph, he::HyperEdge, r1::Root, r2::Root)
    # get (not hyper) vertices of hyperedge
    vs_he_src, vs_he_dst = V(sg, he)

    # find vertices in roots that match source of hyperedge
    src_idx1 = findfirst(v -> v in vs_he_src, vertices(r1))
    src_idx2 = findfirst(v -> v in vs_he_src, vertices(r2))

    # find vertices in roots that match destination of hyperedge (only need to look at vertices neighbouring source idx)
    dst_idx1 = get(vertices(r1), src_idx1 - 1, 0) in vs_he_dst ? src_idx1 - 1 : src_idx1 + 1
    dst_idx2 = get(vertices(r2), src_idx2 - 1, 0) in vs_he_dst ? src_idx2 - 1 : src_idx2 + 1

    # skip switching if hyperedge is at an extremity of either root (switching does nothing)
    if any([idx in [1, length(r1)] for idx in [src_idx1, dst_idx1]]) || any([idx in [1, length(r2)] for idx in [src_idx2, dst_idx2]])
        @debug("`switch!` skipped")
        return nothing
    end

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
