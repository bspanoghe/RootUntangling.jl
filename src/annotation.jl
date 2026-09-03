function write_annotation(filename::String, sg_full::SuperGraph,
        sgs::Vector{<:SuperGraph}, root_systems::Vector{<:Vector{<:Vector{<:Root}}}
    )
    io = open(filename, "w")
    write(io, "Segment_ID,root_ids\n")

    annotation_dict = get_annotation_dict(sgs, root_systems)
    segment_ids = segment_id.(Eₕ₀(sg_full))
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

function append_annotation_dict!(annotation_dict::Dict, sg::SuperGraph, rs::Vector{<:Root})
    for he in Eₕ₀(sg)
        root_idxs = findall(r -> he ∈ Eₕ(r), rs)
        annotation_dict[he] = root_idxs
    end
end

function append_annotation_dict!(annotation_dict::Dict, sg::SuperGraph, root_system::Vector{<:Vector{<:Root}})
    append_annotation_dict!(annotation_dict, sg, reduce(vcat, root_system))
end

function get_annotation_dict(sgs::Vector{<:SuperGraph}, root_systems::Vector{<:Vector{<:Vector{<:Root}}})
    @assert length(sgs) == length(root_systems)
    annotation_dict = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
    num_roots = length.(root_systems)

    for (i, sg) in enumerate(sgs)
        rss = root_systems[i]

        if length(rss) == 1
            rs = rss[1]
            for he in Eₕ₀(sg)
                root_idxs = findall(r -> he ∈ Eₕ(r), rs)
                he_values = get!(annotation_dict, segment_id(he), Tuple{Int64, Int64}[])
                for root_idx in root_idxs
                    new_entry = (sum(num_roots[1:i]), root_idx)
                    new_entry ∈ he_values || (push!(he_values, new_entry))
                end
            end
        else
            for (rs_idx, rs) in enumerate(rss)
                for he in Eₕ₀(sg)
                    root_idxs = findall(r -> he ∈ Eₕ(r), rs)
                    he_values = get!(annotation_dict, segment_id(he), Tuple{Int64, Int64}[])
                    for root_idx in root_idxs
                        new_entry = (sum(num_roots[1:(i-1)]) + rs_idx, root_idx)
                        new_entry ∈ he_values || (push!(he_values, new_entry))
                    end
                end
            end
        end
    end

    return annotation_dict
end

# idk where to put this :(
function get_root_systems(sgs::Vector{SuperGraph{T, U}}, num_roots::Vector{<:Integer};
        optimizer, time_limit, kwargs...
    ) where {T, U}

    root_systems = Vector{Vector{Vector{Root{T, U}}}}(undef, length(sgs))

    for (i, sg) in enumerate(sgs)
        model = solve_rsa(
            sg; optimizer, time_limit, num_roots = num_roots[i], kwargs...
        )

        roots = get_roots(sg, model)
        rss = greedy_switch(sg, model, roots)

        root_systems[i] = rss isa Vector{<:Root} ? [rss] : rss
    end

    return root_systems
end