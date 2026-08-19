function greedy_switch(rg, model, roots; max_tries = 100)
    roots_copy = deepcopy(roots)
    current_tort = tortuosity(roots)
    improving = true

    counter = 0
    while improving
        counter += 1
        counter > max_tries && (@info "Max tries reached"; break)

        overlap_dict = find_overlaps(rg, model, roots_copy)
        switch_dict = get_switch_dict(overlap_dict, roots_copy)
        n = maximum(keys(switch_dict))

        torts = [total_tortuosity(rg, roots_copy, create_switches(n, i), switch_dict) for i in 1:n]
        best_idx = argmin(torts)
        if torts[best_idx] < current_tort
            make_switches!(rg, roots_copy, create_switches(n, best_idx), switch_dict)
        else
            improving = false
        end
    end

    return roots_copy
end

# does a root use one of the edges of a segment
passes_through(r::Root, seg::Vector{<:RootEdge}) = (
    !isdisjoint(srcs(seg), vertices(r)) && !isdisjoint(dsts(seg), vertices(r))
)

# find segments where multiple roots overlap (and can switch)
function find_overlaps(rg::RootGraph, model::JuMP.Model, roots)
    segment_classification_dict = get_segment_classification_dict(rg, model)
    overlap_segments = filter(seg -> imag(segment_classification_dict[seg]) > 1, segments₀(rg)) #! only looks for lateral overlaps

    # map all segments to the roots they are part of
    # discarding roots of length 2 or smaller (switching does nothing)
    overlap_dict = [
        seg => [r for r in roots if passes_through(r, seg) && length(r) > 2]
            for seg in overlap_segments
    ] |> Dict

    # if you have multiple segments in a connected linear path, each with the same amount of roots, only choose one segment from that path
    # reasoning: switching roots on multiple of these consecutive segments has no added effect over switching them on just one
    overlap_segments_subset = [
        get_linear_chains(filter(seg -> length(overlap_dict[seg]) == n, overlap_segments))
            for n in unique(length.(values(overlap_dict)))
    ] |> x -> reduce(vcat, x) .|> first

    overlap_dict_subset = [
        seg => [r for r in roots if passes_through(r, seg)]
            for seg in overlap_segments_subset
    ] |> Dict

    return overlap_dict_subset
end

# group segments into linear chains
function get_linear_chains(segs::Vector{S}) where {S <: Vector{<:RootEdge}}
    chains = Vector{S}[segs[1:1]]
    remaining_segs = segs[2:end]

    while !isempty(remaining_segs)
        chain = chains[end]
        head_idx = findfirst(seg -> !isdisjoint(vertices(chain[1]), vertices(seg)), remaining_segs)
        if !isnothing(head_idx)
            seg = popat!(remaining_segs, head_idx)
            pushfirst!(chain, seg)
            continue
        end

        tail_idx = findfirst(seg -> !isdisjoint(vertices(chain[end]), vertices(seg)), remaining_segs)
        if !isnothing(tail_idx)
            seg = popat!(remaining_segs, tail_idx)
            push!(chain, seg)
            continue
        end

        seg = pop!(remaining_segs)
        push!(chains, [seg])
    end

    return chains
end

function get_switch_dict(overlap_dict, roots)
    switch_pairs = Pair[]
    counter = 0
    for (seg, overlap_roots) in overlap_dict
        for i in eachindex(overlap_roots), j in eachindex(overlap_roots)
            j <= i && continue # order of roots is not important
            counter += 1
            r1_idx = findfirst(r -> r == overlap_roots[i], roots)
            r2_idx = findfirst(r -> r == overlap_roots[j], roots)
            push!(switch_pairs, counter => (seg, r1_idx, r2_idx))
        end
    end

    return Dict(switch_pairs)
end

function total_tortuosity(rg, roots, switches, switch_dict)
    roots_copy = deepcopy(roots)
    make_switches!(rg, roots_copy, switches, switch_dict)
    return tortuosity(roots_copy)
end

create_switches(n::Integer, i::Integer) = [zeros(Bool, i - 1); true; zeros(Bool, n - i)]

function make_switches!(rg, roots, switches, switch_dict)
    for switch in findall(switches)
        seg, r1_idx, r2_idx = switch_dict[switch]
        switch!(rg, seg, roots[r1_idx], roots[r2_idx])
    end

    return nothing
end

# perform a crossing over between two roots
function switch!(rg::RootGraph, seg::Vector{<:RootEdge}, r1::Root, r2::Root)
    # get the vertices on both ends of the segment
    vs_seg_src, vs_seg_dst = srcs(seg), dsts(seg)

    # find vertices in roots that match source of segment
    src_idx1 = findfirst(v -> v in vs_seg_src, vertices(r1))
    src_idx2 = findfirst(v -> v in vs_seg_src, vertices(r2))

    # find vertices in roots that match destination of segment (only need to look at vertices neighbouring source idx)
    dst_idx1 = get(vertices(r1), src_idx1 - 1, 0) in vs_seg_dst ? src_idx1 - 1 : src_idx1 + 1
    dst_idx2 = get(vertices(r2), src_idx2 - 1, 0) in vs_seg_dst ? src_idx2 - 1 : src_idx2 + 1

    # skip switching if segment is at an extremity of either root (switching does nothing)
    if any([idx in [1, length(r1)] for idx in [src_idx1, dst_idx1]]) || any([idx in [1, length(r2)] for idx in [src_idx2, dst_idx2]])
        @debug("`switch!` skipped")
        return nothing
    end

    # check if root orientations match
    orientation_match = (dst_idx1 - src_idx1) == (dst_idx2 - src_idx2)

    # splice roots at segment and switch a half of both
    if orientation_match
        tail1 = splice!(r1.V, dst_idx1:length(r1), r2.V[dst_idx2:end])
        splice!(r2.V, dst_idx2:length(V(r2)), tail1)
    else
        tail1 = splice!(r1.V, dst_idx1:length(r1), reverse(r2.V[1:dst_idx2]))
        splice!(r2.V, 1:dst_idx2, reverse(tail1))
    end

    return nothing
end
