function write_annotation(filename::String, rg_full::RootGraph, annotation_dict::Dict)
    io = open(filename, "w")
    write(io, "Segment_ID,root_ids,real_length\n")

    segment_ids = segment_id.(E₀(rg_full))
    sort!(segment_ids)

    for segment_id in segment_ids
        root_ids = get(annotation_dict, segment_id, Int64[])
        root_id_string = isempty(root_ids) ? "-" :
            *([get_string_id(id) * ";" for id in root_ids]...)[1:end-1] # remove trailing `,`
        segment_id_string = get_string_id(segment_id)
        write(io, segment_id_string * "," * root_id_string * ",\n")
    end
    close(io)
end

get_string_id(id::Integer) = string(id)
get_string_id(id::Tuple{Int64, Int64}) = string(id[1]) * "." * string(id[2])

# get `edge segment id => [roots containing edge]`, e.g. `1 => [(1, 1)], 2 => [(1, 1), (1, 2)]`
function get_annotation_dict(rgs::Vector{<:RootGraph}, root_clusters::Vector{<:Vector{<:RootSystem}})
    @assert length(rgs) == length(root_clusters)
    annotation_dict = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
    num_roots = length.(root_clusters)

    for (i, rg) in enumerate(rgs)
        rss = root_clusters[i]

        if length(rss) == 1
            rs = rss[1]
            for re in E₀(rg)
                root_idxs = findall(r -> re ∈ edges(r), roots(rs))
                re_values = get!(annotation_dict, segment_id(re), Tuple{Int64, Int64}[])
                for root_idx in root_idxs
                    new_entry = (sum(num_roots[1:i]), root_idx)
                    new_entry ∈ re_values || (push!(re_values, new_entry))
                end
            end
        else
            for (rs_idx, rs) in enumerate(rss)
                for re in E₀(rg)
                    root_idxs = findall(r -> re ∈ edges(r), roots(rs))
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
function get_root_clusters(rgs::Vector{RootGraph{T, U}}, pg::PreGraph{T, U}, num_roots::Vector{<:Integer};
        optimizer, time_limit, kwargs...
    ) where {T, U}

    root_clusters = Vector{Vector{RootSystem{T, U}}}(undef, length(rgs))

    for (i, rg) in enumerate(rgs)
        model = solve_rsa(
            rg; optimizer, time_limit, num_roots = num_roots[i], kwargs...
        )

        roots = get_rootsystems(rg, pg, model)
        rss = greedy_switch(roots)
        root_clusters[i] = rss
    end

    return root_clusters
end