function write_annotation(filename::String, rg_full::RootGraph,
        rgs::Vector{<:RootGraph}, root_systems::Vector{<:Vector{<:Vector{<:Root}}}
    )
    io = open(filename, "w")
    write(io, "Segment_ID,root_ids\n")

    annotation_dict = get_annotation_dict(rgs, root_systems)
    segment_ids = segment_id.(Eₕ₀(rg_full))
    sort!(segment_ids)

    for segment_id in segment_ids
        root_ids = get(annotation_dict, segment_id, Int64[])
        root_id_string = isempty(root_ids) ? "-" :
            *([get_string_id(id) * ";" for id in root_ids]...)[1:end-1] # remove trailing `,`
        segment_id_string = get_string_id(segment_id)
        write(io, segment_id_string * "," * root_id_string * "\n")
    end
    close(io)
end

get_string_id(id::Integer) = string(id)
get_string_id(id::Tuple{Int64, Int64}) = string(id[1]) * "." * string(id[2])

function append_annotation_dict!(annotation_dict::Dict, rg::RootGraph, rs::Vector{<:Root})
    for re in Eₕ₀(rg)
        root_idxs = findall(r -> re ∈ Eₕ(r), rs)
        annotation_dict[re] = root_idxs
    end
end

function append_annotation_dict!(annotation_dict::Dict, rg::RootGraph, root_system::Vector{<:Vector{<:Root}})
    append_annotation_dict!(annotation_dict, rg, reduce(vcat, root_system))
end

function get_annotation_dict(rgs::Vector{<:RootGraph}, root_systems::Vector{<:Vector{<:Vector{<:Root}}})
    @assert length(rgs) == length(root_systems)
    annotation_dict = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
    num_roots = length.(root_systems)

    for (i, rg) in enumerate(rgs)
        rss = root_systems[i]

        if length(rss) == 1
            rs = rss[1]
            for re in Eₕ₀(rg)
                root_idxs = findall(r -> re ∈ Eₕ(r), rs)
                re_values = get!(annotation_dict, segment_id(re), Tuple{Int64, Int64}[])
                for root_idx in root_idxs
                    new_entry = (sum(num_roots[1:i]), root_idx)
                    new_entry ∈ re_values || (push!(re_values, new_entry))
                end
            end
        else
            for (rs_idx, rs) in enumerate(rss)
                for re in Eₕ₀(rg)
                    root_idxs = findall(r -> re ∈ Eₕ(r), rs)
                    re_values = get!(annotation_dict, segment_id(re), Tuple{Int64, Int64}[])
                    for root_idx in root_idxs
                        new_entry = (sum(num_roots[1:(i-1)]) + rs_idx, root_idx)
                        new_entry ∈ re_values || (push!(re_values, new_entry))
                    end
                end
            end
        end
    end

    return annotation_dict
end

# idk where to put this :(
function get_root_systems(rgs::Vector{RootGraph{T, U}}, num_roots::Vector{<:Integer};
        optimizer, time_limit, kwargs...
    ) where {T, U}

    root_systems = Vector{Vector{Vector{Root{T, U}}}}(undef, length(rgs))

    for (i, rg) in enumerate(rgs)
        model = solve_rsa(
            rg; optimizer, time_limit, num_roots = num_roots[i], kwargs...
        )

        roots = get_roots(rg, model)
        rss = greedy_switch(rg, model, roots)

        root_systems[i] = rss isa Vector{<:Root} ? [rss] : rss
    end

    return root_systems
end