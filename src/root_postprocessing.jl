function greedy_switch(rg::RootGraph, model::JuMP.Model, roots::Vector{<:Root};
        max_tries = 100, f_obj::Function = roughness
    )
    if isempty(find_overlaps(rg, model, roots))
        @info "No overlaps found"
        return roots
    end

    roots_copy = deepcopy(roots)
    current_f = f_obj(roots)
    improving = true

    counter = 0
    while improving
        counter += 1
        counter > max_tries && (@info "Max tries reached"; break)

        overlap_dict = find_overlaps(rg, model, roots_copy)
        switch_dict = get_switch_dict(overlap_dict, roots_copy)      
        n = maximum(keys(switch_dict))

        fs = [evaluate_objective(rg, roots_copy, f_obj, create_switches(n, i), switch_dict) for i in 1:n]
        best_idx = argmin(fs)
        if fs[best_idx] < current_f
            make_switches!(rg, roots_copy, create_switches(n, best_idx), switch_dict)
            current_f = fs[best_idx]
        else
            improving = false
        end
    end

    return roots_copy
end

function greedy_switch(rg::RootGraph, model::JuMP.Model, root_systems::Vector{<:Vector{<:Root}};
        max_tries = 100, f_obj::Function = roughness
    )
    
    roots_tangled = greedy_switch(rg, model, reduce(vcat, root_systems); max_tries, f_obj)
    return separate_root_systems(rg, get_re_classification_dict(rg, model), roots_tangled)
end

# does a root use one of the edges of a segment
are_overlapping(rg::RootGraph, re::RootEdge, r::Root) = (
    !all([isdisjoint(vs, vertices(r)) for vs in vertices(Vₕ(rg, re))])
)

# find hyperedges where multiple lateral roots overlap (and can switch)
# primary root is not allowed to swap because it can break the assumption that only the primary root can split
function find_overlaps(rg::RootGraph, model::JuMP.Model, roots)
    re_classification_dict = get_re_classification_dict(rg, model)
    overlap_res = filter(re -> imag(re_classification_dict[re]) > 1, Eₕ₀(rg))

    # map all hyperedges to the roots they are part of
    # discarding roots of length 2 or smaller (switching does nothing)
    overlap_dict = [
        re => [r for r in roots if are_overlapping(rg, re, r) && length(r) > 2]
            for re in overlap_res
    ] |> Dict{RootEdge, Vector{<:Root}}
    
    isempty(overlap_dict) && (return Dict{RootEdge, Vector{<:Root}}())

    # if you have multiple edges in a connected linear path, each with the same amount of roots, only choose one edge from that path
    # reasoning: switching roots on multiple of these consecutive edges has no added effect over switching them on just one
    overlap_res_subset = [
        get_linear_chains(filter(re -> length(overlap_dict[re]) == n, overlap_res))
            for n in unique(length.(values(overlap_dict)))
    ] |> x -> reduce(vcat, x) .|> first

    overlap_dict_subset = [
        re => [r for r in roots if are_overlapping(rg, re, r)]
            for re in overlap_res_subset
    ] |> Dict{RootEdge, Vector{<:Root}}

    return overlap_dict_subset
end

# group edges into linear chains
function get_linear_chains(res::Vector{<:RootEdge})
    chains = Vector{RootEdge}[res[1:1]]
    remaining_res = res[2:end]

    while !isempty(remaining_res)
        chain = chains[end]
        head_idx = findfirst(re -> !isdisjoint(vertices(chain[1]), vertices(re)), remaining_res)
        if !isnothing(head_idx)
            re = popat!(remaining_res, head_idx)
            pushfirst!(chain, re)
            continue
        end

        tail_idx = findfirst(re -> !isdisjoint(vertices(chain[end]), vertices(re)), remaining_res)
        if !isnothing(tail_idx)
            re = popat!(remaining_res, tail_idx)
            push!(chain, re)
            continue
        end

        re = pop!(remaining_res)
        push!(chains, [re])
    end

    return chains
end

function get_switch_dict(overlap_dict, roots)
    switch_pairs = Pair[]
    counter = 0
    for (re, overlap_roots) in overlap_dict
        for i in eachindex(overlap_roots), j in eachindex(overlap_roots)
            j <= i && continue # order of roots is not important
            counter += 1
            r1_idx = findfirst(r -> r == overlap_roots[i], roots)
            r2_idx = findfirst(r -> r == overlap_roots[j], roots)
            push!(switch_pairs, counter => (re, r1_idx, r2_idx))
        end
    end

    return Dict(switch_pairs)
end

function evaluate_objective(rg::RootGraph, roots::Vector{<:Root}, f_obj::Function, switches::Vector{Bool}, switch_dict::Dict)
    roots_copy = deepcopy(roots)
    make_switches!(rg, roots_copy, switches, switch_dict)
    return f_obj(roots_copy)
end

create_switches(n::Integer, i::Integer) = [zeros(Bool, i - 1); true; zeros(Bool, n - i)]

function make_switches!(rg, roots, switches, switch_dict)
    for switch in findall(switches)
        re, r1_idx, r2_idx = switch_dict[switch]
        switch!(rg, re, roots[r1_idx], roots[r2_idx])
    end

    return nothing
end

# perform a crossing over between two roots
function switch!(rg::RootGraph, re::RootEdge, r1::Root, r2::Root)
    # get (not hyper) vertices of hyperedge
    vs_re_src, vs_re_dst = vertices.(Vₕ(rg, re)) # `.` to get separately for both hypervertices

    # find vertices in roots that match source of hyperedge
    src_idx1 = findfirst(v -> v in vs_re_src, vertices(r1))
    src_idx2 = findfirst(v -> v in vs_re_src, vertices(r2))

    # find vertices in roots that match destination of hyperedge (only need to look at vertices neighbouring source idx)
    dst_idx1 = get(vertices(r1), src_idx1 - 1, 0) in vs_re_dst ? src_idx1 - 1 : src_idx1 + 1
    dst_idx2 = get(vertices(r2), src_idx2 - 1, 0) in vs_re_dst ? src_idx2 - 1 : src_idx2 + 1

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
