struct ValRoot{T, U}
    is_primary::Bool
    edges::Vector{RootEdge{T, U}}
    xs::Vector{U}
    ys::Vector{U}
end
is_primary(vr::ValRoot) = vr.is_primary
edges(vr::ValRoot) = vr.edges
xs(vr::ValRoot) = vr.xs
ys(vr::ValRoot) = vr.ys

Base.length(vr::ValRoot) = length(edges(vr))
Makie.convert_arguments(::Type{<:Scatter}, vr::ValRoot) = (xs(vr), ys(vr))

curve_length(vr::ValRoot) = curve_length(xs(vr), ys(vr))

struct ValRootSystem{T, U}
    primary::ValRoot{T, U}
    laterals::Vector{ValRoot{T, U}}
end
primary(vrs::ValRootSystem) = vrs.primary
laterals(vrs::ValRootSystem) = vrs.laterals
roots(vrs::ValRootSystem) = [primary(vrs); laterals(vrs)]
Base.length(vrs::ValRootSystem) = length(roots(vrs))
xs(vrs::ValRootSystem) = reduce(vcat, xs.(roots(vrs)))
ys(vrs::ValRootSystem) = reduce(vcat, ys.(roots(vrs)))

function get_validation_data(filename_GT::String, filename_segments::String, filename_vertices::String;
        cluster_threshold, loop_threshold, reverse_y, angle_correction
    )
    
    # get original (not cleaned) segments
    vertex_data_dict, edge_data_dict = get_data_dicts(filename_segments, filename_vertices; reverse_y, angle_correction)
    metavertices = getmetavertices(vertex_data_dict)
    segments = getsegments(metavertices, vertex_data_dict, edge_data_dict)
    sort_segment_coordinates!(segments, metavertices)

    # get how segments map to roots in the ground truth
    valdata = read_data(filename_GT, :Segment_ID, missing_string = "-", parse_strings = false) |>
        x -> [parse(Int64, key) => x[key] for key in keys(x)] |> Dict
    root_segment_mappings = get_root_segment_mappings(valdata)

    rg = get_rootgraph(filename_segments, filename_vertices; cluster_threshold, loop_threshold, reverse_y, angle_correction)

    vrss = get_valrootsystem.([rg], [segments], root_segment_mappings)

    root_segmentation_errors = [
        filter(!ismissing, Any[
            valdata[seg_id][:real_length]
            for seg_id in reduce(vcat, values(root_segment_mapping))
        ]) .|> (x -> parse(Int64, x)) |> sum
        for root_segment_mapping in root_segment_mappings
    ]

    return vrss, root_segmentation_errors
end

function get_root_segment_mappings(valdata::Dict)
    id_mapping = [
        key => [
            parse.(Int64, root_id) 
            for root_id in split.(valdata[key][:root_ids], '.')
        ] |>
            x -> x isa Vector{Int64} ? [(x...,)] : [(_x...,) for _x in x]
        for key in keys(valdata)
        if !ismissing(valdata[key][:root_ids])
    ] |> Dict

    root_ids = values(id_mapping) |>
        x -> reduce(vcat, x) |>
        unique

    root_mapping = [
        id => [key for key in keys(id_mapping) if id ∈ id_mapping[key]]
        for id in root_ids
    ]

    rs_ids = unique(first.(first.(root_mapping))) |> sort
    root_segment_mappings = [
        [
            mapping[1][2] => mapping[2]
            for mapping in root_mapping
            if mapping[1][1] == rs_id
        ] |> Dict
        for rs_id in rs_ids
    ]
    @assert all(allunique(keys(root_segment_mappings[i])) for i in eachindex(root_segment_mappings))

    return root_segment_mappings
end

function get_valrootsystem(rg::RootGraph{T, U}, segments, root_segment_mapping::Dict) where {T, U}
    primary = ValRoot(
        true,
        filter(e -> e.segment_id ∈ root_segment_mapping[1], E₀(rg)),
        reduce(vcat, xs.(filter(s -> s.id ∈ root_segment_mapping[1], segments)), init = U[]),
        reduce(vcat, ys.(filter(s -> s.id ∈ root_segment_mapping[1], segments)), init = U[]),
    )
    laterals = [
        ValRoot(
            false,
            filter(e -> e.segment_id ∈ root_segment_mapping[i], E₀(rg)),
            reduce(vcat, xs.(filter(s -> s.id ∈ root_segment_mapping[i], segments)), init = U[]),
            reduce(vcat, ys.(filter(s -> s.id ∈ root_segment_mapping[i], segments)), init = U[]),
        )
        for i in filter(x -> x != 1, keys(root_segment_mapping))
    ]

    sort_root_coordinates!(primary, rg)
    for lateral in laterals
        sort_root_coordinates!(lateral, rg)
    end

    return ValRootSystem(primary, laterals)
end

function sort_root_coordinates!(vr::ValRoot, rg::RootGraph)
    vs = collect(Iterators.flatten(vertices.(edges(vr))))
    isempty(vs) && return nothing
    v_start = vs[findfirst(v -> !(v in duplicate_elements(vs)), vs)]
    sort_coordinates!(vr.xs, vr.ys, V(rg, v_start).x, V(rg, v_start).y)

    return nothing
end

function fix_groundtruth(file::String)
    lines = readlines(file)
    filter!(!isempty, lines)
    isempty(split(lines[1], ',')[end]) && (lines[1] *= "real_length")

    datamat = [
        split(line, ',') |> 
            x -> length(x) == 2 ? [x; ""] : x
        for line in lines
    ] |> x -> reduce(hcat, x) |> permutedims

    # (potentially) fix duplicate ids
    ids = datamat[2:end, 1]
    if !allunique(ids)
        @warn("Duplicate ids found in ground truth. Attempting repair.")
        duplicate_ids = unique(duplicate_elements(ids))
        delete_idxs = Int64[]
        for id in duplicate_ids
            row_idxs = findall(datamat[:, 1] .== id)
            rows = [datamat[r, :] for r in row_idxs]
            if allequal(rows)
                @debug "id $(id): All rows equal"
                append!(delete_idxs, row_idxs[2:end])
            else
                @debug "id $(id): Unequal rows\n$(rows)"
                truerow_idx = findfirst(
                    r -> all([
                        all([occursin(row[i], r[i]) for i in eachindex(r)])
                        for row in rows
                    ]), 
                    rows
                )
                isnothing(truerow_idx) && @error "id $(id): Unequal rows\n$(rows)\nTrue data is ambiguous"
                @debug "True row found\n$(rows[truerow_idx])"
                append!(delete_idxs, row_idxs[[1:truerow_idx-1; truerow_idx+1:end]])
            end
        end
        datamat = datamat[.!in.(1:size(datamat, 1), [delete_idxs]), :]
    end

    # correct root annotations to have sequential root idxs
    ann_correction_dict = Dict{Regex, String}()
    root_anns = reduce(vcat, split.(datamat[2:end, 2], ';')) |>
        x -> filter(_x -> _x != "-", x) |>
        x -> split.(x, '.') .|>
        x -> isempty(x) ? x : parse.(Int64, x)

    for rs_idx in unique(first.(root_anns))
        current_root_anns = filter(x -> x[1] == rs_idx, root_anns) .|> last |> unique |> sort
        for (i, og_idx) in enumerate(current_root_anns)
            i != og_idx && (
                ann_correction_dict[Regex("$(rs_idx).$(og_idx)(?![0-9])")] = "$(rs_idx).$(i)"
            )
        end
    end

    isempty(ann_correction_dict) || @info "$(length(ann_correction_dict)) root indices corrected to ensure sequential root indices"
    for i in 2:size(datamat, 1)
        datamat[i, 2] = replace(datamat[i, 2], collect(ann_correction_dict)...)
    end

    # convert back to text and write
    newlines = [*(r .* ','...)[1:end-1] for r in eachrow(datamat)]
    newtext = *((newlines .* '\n')...)[1:end-1]

    write(file, newtext)

    return nothing
end