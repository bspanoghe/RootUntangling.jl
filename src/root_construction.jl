"""
    get_rootsystems(rg::RootGraph, model::JuMP.Model)

Extract the rootsystems from a graph `rg` and its solution contained in `model`.
"""
function get_rootsystems(rg::RootGraph{T, U}, model::JuMP.Model) where {T, U}
    
    # get primary roots
    c_counts = Dict(E₂(rg) .=> round.(Int64, value.(model[:cp])))
    e_counts = Dict(E(rg) .=> round.(Int64, value.(model[:ep])))
    e₊_counts = Dict(E(rg) .=> round.(Int64, value.(model[:ep₊])))
    e₋_counts = Dict(E(rg) .=> round.(Int64, value.(model[:ep]) - value.(model[:ep₊])))

    @debug "It begins"
    primary_fragments = fragment(rg, true, c_counts, e_counts, e₊_counts, e₋_counts)
    @debug "The first gate"
    stitch_fragments!(primary_fragments)
    @debug "The second gate"
    primary_roots = get_roots(primary_fragments)
    @debug "The third gate"

    # get lateral roots
    c_counts = Dict(E₂(rg) .=> round.(Int64, value.(model[:cl])))
    e_counts = Dict(E(rg) .=> round.(Int64, value.(model[:el])))
    e₊_counts = Dict(E(rg) .=> round.(Int64, value.(model[:el₊])))
    e₋_counts = Dict(E(rg) .=> round.(Int64, value.(model[:el]) - value.(model[:el₊])))

    lateral_fragments = fragment(rg, true, c_counts, e_counts, e₊_counts, e₋_counts)
    @debug "The fourth gate"
    stitch_fragments!(lateral_fragments)
    @debug "The fifth gate"
    lateral_roots = get_roots(lateral_fragments)
    @debug "The final gate"

    return separate_root_systems(rg, primary_roots, lateral_roots)
end

function fragment(rg::RootGraph{T, U}, is_primary::Bool,
        c_counts::Dict, e_counts::Dict, e₊_counts::Dict, e₋_counts::Dict
    ) where {T, U}

    fragments = RootFragment[]

    for c in filter(c -> c_counts[c] > 0, E₂(rg))
        v_shared = shared_vertex(c)

        if c_counts[c] == e_counts[c[1]] # first edge fully explains connection

            e = c[1]
            for i in 1:e₊_counts[e] # root follows direction of edge
                if dst(e) == v_shared # does edge go toward shared vertex?
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = false)
                else
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = true)
                end                    
            end

            for i in 1:e₋_counts[e] # root goes against direction of edge
                if dst(e) != v_shared # flip direction
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = false)
                else
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = true)
                end        
            end

        elseif c_counts[c] == e_counts[c[2]] # second edge fully explains connection

            e = c[2]
            for i in 1:e₊_counts[e] # root follows direction of edge
                if src(e) == v_shared # does edge go away from shared vertex?
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = false)
                else
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = true)
                end
            end

            for i in 1:e₋_counts[e] # root goes against direction of edge
                if src(e) != v_shared # flip direction
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = false)
                else
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = true)
                end
            end

        else # neither edge fully explains connection

            @warn "OHHHH NOOOOOOOOOOOOOO"
            for i in 1:c_counts[c]
                push!(fragments, UndirectedRootFragment(is_primary, c))
            end

        end
    end

    return fragments
end

shared_vertex(c::Vector{<:RootEdge}) = vertices(c[1])[findfirst(v -> in(v, vertices(c[2])), vertices(c[1]))]
function add_directed_fragment!(fragments::Vector{<:RootFragment}, is_primary::Bool, c::Vector{<:RootEdge}, v_shared; switch_edges::Bool)
    if !switch_edges
        push!(fragments, DirectedRootFragment(
            is_primary, 
            [
                RootArc(c[1], keep_order = (dst(c[1]) == v_shared)),
                RootArc(c[2], keep_order = (src(c[2]) == v_shared))
            ])
        )
    else
        push!(fragments, DirectedRootFragment(is_primary,
            [
                RootArc(c[2], keep_order = (dst(c[2]) == v_shared)),
                RootArc(c[1], keep_order = (src(c[1]) == v_shared))
            ])
        )
    end

    return nothing
end

function stitch_fragments!(fragments::Vector{<:RootFragment})
    growing = true
    # keep going until all root fragments cant grow anymore
    while growing
        growing = false
        for (i, fragment) in enumerate(fragments)
            
            f_end_idxs = findall(f -> connects_to_end(fragment, f), fragments)
            if are_options_unambiguous(fragments[f_end_idxs]) && fragment != fragments[f_end_idxs[1]] 
                stitch_to_end!(fragment, fragments[f_end_idxs[1]])
                deleteat!(fragments, f_end_idxs[1])
                
                growing = true
            end

            f_start_idxs = findall(f -> connects_to_start(fragment, f), fragments)
            if are_options_unambiguous(fragments[f_start_idxs]) && fragment != fragments[f_start_idxs[1]]
                stitch_to_start!(fragment, fragments[f_start_idxs[1]]) # connect one of them
                deleteat!(fragments, f_start_idxs[1])

                growing = true
            end

        end
    end

    return nothing
end

are_options_unambiguous(xs) = (
    (length(xs) == 1) || # there is only one option
        (length(xs) > 1 && allequal(xs)) # all options are equal
)

connects_to_end(f1::DirectedRootFragment, f2::DirectedRootFragment) = edges(f1)[end] == edges(f2)[1]
connects_to_start(f1::DirectedRootFragment, f2::DirectedRootFragment) = edges(f1)[1] == edges(f2)[end]

stitch_to_end!(f1::DirectedRootFragment, f2::DirectedRootFragment) = append!(f1.edges, f2.edges[2:end])
stitch_to_start!(f1::DirectedRootFragment, f2::DirectedRootFragment) = prepend!(f1.edges, f2.edges[1:end-1])

function get_roots(fs::Vector{<:RootFragment})
    roots = Root[]
    # roots consisting of a single root fragment
    complete_fragments = filter(is_fullgrown, fs)
    append!(roots, Root[SimpleRoot(is_primary(f), f) for f in complete_fragments])

    # roots consisting of multiple root fragments
    incomplete_fragments = filter(!is_fullgrown, fs)

    # keep going until all incomplete root fragments are used
    while !isempty(incomplete_fragments)

        # instantiate new root
        fragment = incomplete_fragments[1]
        deleteat!(incomplete_fragments, 1)
        root = CompositeRoot(is_primary(fragment), [fragment])

        growing = true # needed in case of loops :[
        # grow until complete
        while !is_fullgrown(root) && growing
            fragment = fragments(root)[end]
            f_end_idxs = findall(f -> connects_to_end(fragment, f), incomplete_fragments)
            if !isempty(f_end_idxs)
                push!(fragments(root), incomplete_fragments[f_end_idxs[1]])
                deleteat!(incomplete_fragments, f_end_idxs[1])
            end

            fragment = fragments(root)[1]
            f_start_idxs = findall(f -> connects_to_start(fragment, f), incomplete_fragments)
            if !isempty(f_start_idxs)
                pushfirst!(fragments(root), incomplete_fragments[f_start_idxs[1]])
                deleteat!(incomplete_fragments, f_start_idxs[1])
            end

            growing = !(isempty(f_end_idxs) && isempty(f_start_idxs))
        end
        push!(roots, root)
    end

    return roots
end

# divide roots into separate root systems
function separate_root_systems(rg::RootGraph, primary_roots::Vector{<:Root}, lateral_roots::Vector{<:Root})
    lateral_match_dict = Dict{Root, Int64}()

    for lateral_root in lateral_roots
        start_edge = edges(lateral_root)[1]
        if vertices(start_edge)[1] == -3 # lateral splits from a primary root
            lateral_match_dict[lateral_root] = findfirst(
                r -> vertices(start_edge)[2] in vertices(r), # first vertex is -3
                primary_roots
            )
        else # lateral appears straight up out of nowhere
            lateral_match_dict[lateral_root] = findmin(
                pr -> distance(rg, V₀(rg, lateral_root)[1], pr), 
                primary_roots
            )[2] # findmin returns (element, idx)
        end
    end

    return [
        RootSystem(primary_root, filter(r -> lateral_match_dict[r] == i, lateral_roots)) 
        for (i, primary_root) in enumerate(primary_roots)
    ]
end

function sort_root_system!(rs::Vector{<:Root})
    sort!(rs, by = is_primary, rev = true)
    if length(rs) > 1
        @assert !is_primary(rs[2]) "Root systems should only contain one primary root"
        sort!(@view(rs[2:end]), by = curve_length, rev = true)
    end

    return nothing
end