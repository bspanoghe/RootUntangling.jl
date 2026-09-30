function greedy_switch(rg::RootGraph, rs::RootSystem; max_tries = 100, f_obj::Function = weighted_tortuosity)

    composite_laterals = filter(r -> r isa CompositeRoot, laterals(rs))
    if isempty(composite_laterals)
        @info "Root system has no overlaps"
        return rs
    end
    composite_laterals_switched = greedy_switch(rg, composite_laterals; max_tries, f_obj)
    simple_laterals = filter(r -> r isa SimpleRoot, laterals(rs))

    return RootSystem(
        primary(rs),
        sort([simple_laterals; composite_laterals_switched], by = r -> curve_length(rg, r), rev = true)
    )
end

function greedy_switch(rg::RootGraph, rss::Vector{<:RootSystem}; max_tries = 100, f_obj::Function = weighted_tortuosity)
    
    composite_laterals = reduce(vcat, [filter(r -> r isa CompositeRoot, laterals(rs)) for rs in rss])
    if isempty(composite_laterals)
        @info "Root systems have no overlaps"
        return rss
    end
    composite_laterals_switched = greedy_switch(rg, composite_laterals; max_tries, f_obj)
    simple_laterals = reduce(vcat, [filter(r -> r isa SimpleRoot, laterals(rs)) for rs in rss])

    return separate_root_systems(rg, primary.(rss), [simple_laterals; composite_laterals_switched])
end

function greedy_switch(rg::RootGraph, roots::Vector{<:Root}; max_tries::Integer, f_obj::Function)
    roots_copy = deepcopy(roots)
    current_f = f_obj(rg, roots)
    improving = true

    counter = 0
    while improving
        counter += 1
        @debug(counter)
        counter > max_tries && (@info "Max tries reached"; break)

        overlaps = find_overlaps(roots_copy)

        fs = [evaluate_objective(rg, roots_copy, f_obj, overlap) for overlap in overlaps]
        best_idx = argmin(fs)
        if fs[best_idx] < current_f
            switch!(roots_copy; overlaps[best_idx]...)
            current_f = fs[best_idx]
        else
            improving = false
        end
    end

    return roots_copy
end

# find what roots overlap and at which fragments
function find_overlaps(roots::Vector{<:Root})
    overlaps = [
        (r1 = r1, r2 = r2, f1 = f1, f2 = f2)
        for r1 in eachindex(roots) for r2 in eachindex(roots)[(r1+1):end]
        for f1 in eachindex(fragments(roots[r1])) for f2 in eachindex(fragments(roots[r2]))
        if are_overlapping(fragments(roots[r1])[f1], fragments(roots[r2])[f2])
    ]
    
    return overlaps
end

# do two root fragments overlap
are_overlapping(rf1::RootFragment, rf2::RootFragment) = (edges(rf1)[1] == edges(rf2)[1]) ||
    (edges(rf1)[end] == edges(rf2)[end])

function evaluate_objective(rg::RootGraph, roots::Vector{<:Root}, f_obj::Function, overlap::NamedTuple)
    roots_copy = deepcopy(roots)
    switch!(roots_copy; overlap...)

    return f_obj(rg, roots_copy)
end

# perform a crossing over between two roots at specified root fragments
function switch!(rs::Vector{<:Root}; r1::T, r2::T, f1::T, f2::T) where {T <: Integer}
    tail1 = splice!(rs[r1].fragments, f1:length(rs[r1].fragments), fragments(rs[r2])[f2:end])
    splice!(rs[r2].fragments, f2:length(rs[r2].fragments), tail1)

    return nothing
end