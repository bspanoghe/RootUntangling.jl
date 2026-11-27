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

function get_valrootsystems(filename_GT::String, pg::PreGraph, rg::RootGraph)
    # read ground truth `segment_id` => `[(rootsys_idx, rootnum_idx), (...), ...]`
    annotation_dict = read_annotation_dict(filename_GT)

    # cluster segment ids by rootsystem and root
    ## invert annotation dict `(rootsys_idx, rootnum_idx) => segment_ids`
	inverted_annotation_dict = Dict{valtype(annotation_dict).parameters[1], Vector{keytype(annotation_dict)}}()
    for (k, v) in annotation_dict
        for segment_id in v
            push!(
                get!(inverted_annotation_dict, segment_id, Int64[]),
                k
            )
        end
    end
    ## expand dict to `rootsys_idx` => `rootnum_idx` => segment_ids
    rs_idxs = unique(first.(keys(inverted_annotation_dict)))
    root_seg_dicts = [
        Dict([
            rootnum_idx => inverted_annotation_dict[(rs_idx, rootnum_idx)]
            for rootnum_idx in unique(last.(filter(k -> k[1] == rs_idx, keys(inverted_annotation_dict))))
        ])
        for rs_idx in rs_idxs
    ]

    # make validation rootsystems from clustered rootedges
    return [get_valrootsystem(root_seg_dict, pg, rg) for root_seg_dict in root_seg_dicts]
end

function get_valrootsystem(root_seg_dict::Dict, pg::PreGraph{T, U}, rg::RootGraph{T, U}) where {T, U}
    primary = ValRoot(
        true,
        filter(re -> og_segment_id(re) ∈ root_seg_dict[1], E₀(rg)),
        reduce(vcat, xs.(filter(s -> og_id(s) ∈ root_seg_dict[1], pg.segments)), init = U[]),
        reduce(vcat, ys.(filter(s -> og_id(s) ∈ root_seg_dict[1], pg.segments)), init = U[])
    )
    laterals = [
        ValRoot(
            false,
            filter(re -> og_segment_id(re) ∈ root_seg_dict[i], E₀(rg)),
            reduce(vcat, xs.(filter(s -> og_id(s) ∈ root_seg_dict[i], pg.segments)), init = U[]),
            reduce(vcat, ys.(filter(s -> og_id(s) ∈ root_seg_dict[i], pg.segments)), init = U[])
        )
        for i in filter(x -> x != 1, keys(root_seg_dict))
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