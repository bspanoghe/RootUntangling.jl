function fix_groundtruth(file::String)
    lines = readlines(file)
    filter!(!isempty, lines)
    newlines = [
        length(split(line, ',')) == 2 ? line * ',' : line
        for line in lines
    ]
    isempty(split(newlines[1], ',')[end]) && (newlines[1] *= "real_length")
    newtext = *((newlines .* '\n')...)[1:end-1]

    write(file, newtext)

    return nothing
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

function rootplot!(ax::Makie.Axis, vrs::ValRootSystem; kwargs...)
    r = primary(vrs)
    scatter!(ax, r; color = HSV(0, 1, 0), label = "0", markersize = 2, kwargs...)

    for (i, r) in enumerate(laterals(vrs))
        color = HSV(range(0, 360, length = length(vrs))[i], 1, 0.75)
        scatter!(ax, r; color, label = "$i", markersize = 2, kwargs...)
    end
end

function rootplot!(ax::Makie.Axis, vrss::Vector{<:ValRootSystem}; kwargs...)
    for (i, vrs) in enumerate(vrss)
        r = primary(vrs)
        scatter!(ax, r; color = HSV(0, 1, 0), label = "$(i)-0", markersize = 2, kwargs...)

        for (j, r) in enumerate(laterals(vrs))
            color = HSV(range(0, 360, length = length(vrs))[j], 1, 0.75)
            scatter!(ax, r; color, label = "$(i)-$(j)", markersize = 2, kwargs...)
        end
    end
end

function rootplot(vrs::ValRootSystem; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, vrs; line_kwargs...)

    return f
end

function rootplot(vrss::Vector{<:ValRootSystem}; size = (1200, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, vrss; line_kwargs...)

    return f
end

function rootplot(vrs::ValRootSystem, rs::RootSystem; size = (1200, 600), kwargs...)
    f = Figure(; size)
    rootplot!(Axis(f[1, 1], aspect = DataAspect()), vrs)
    f[1, 2] = Legend(f, f.content[end])
    rootplot!(Axis(f[1, 3], aspect = DataAspect()), rs)
    f[1, 4] = Legend(f, f.content[end])
    return f
end

add_dummy_roots!(rs::RootSystem, vrs::ValRootSystem) = add_dummy_roots!(laterals(rs), laterals(vrs))

function add_dummy_roots!(rs::Vector{<:Root}, vrs::Vector{<:ValRoot})
    n_diff = length(rs) - length(vrs)
    n_diff < 0 && add_dummy_roots!(rs, abs(n_diff))
    n_diff > 0 && add_dummy_roots!(vrs, n_diff)

    return nothing
end

function add_dummy_roots!(roots::Vector{Union{CompositeRoot{T, U}, SimpleRoot{T, U}}}, n_diff::Integer) where {T, U}
    append!(roots, fill(SimpleRoot(false, RootFragment(false, RootArc{T, U}[])), n_diff))

    return nothing
end

function add_dummy_roots!(roots::Vector{ValRoot{T, U}}, n_diff::Integer) where {T, U}
    append!(roots, fill(ValRoot(false, RootEdge{T, U}[], U[], U[]), n_diff))

    return nothing
end

function get_loss(rs::RootSystem, vrs::ValRootSystem, root_mapping::Vector{<:Integer}, metric::Function;
        verbose::Bool = false
    )
    #! prevent deepcopy here for efficiency
    rs, vrs = deepcopy(rs), deepcopy(vrs)
    add_dummy_roots!(rs, vrs)

    # primary root
    primary_dist = get_loss(primary(rs), primary(vrs), metric)

    # laterals
    lateral_dist = 0
    for i in eachindex(root_mapping)
        lateral_dist += get_loss(laterals(rs), laterals(vrs), i, root_mapping[i], metric)
        verbose && @info "Loss in root $(i): $(get_loss(laterals(rs), laterals(vrs), i, root_mapping[i], metric))"
    end

    return primary_dist + lateral_dist
end

get_loss(rs::Vector{<:Root}, rvs::Vector{<:ValRoot}, i::Integer, j::Integer, metric::Function) = get_loss(rs[i], rvs[j], metric)

function get_loss(r::Root, rv::ValRoot, metric::Function)
    return length(full_setdiff(metric(r), metric(rv)))
end

full_setdiff(x1::AbstractArray, x2::AbstractArray) = [setdiff(x1, x2); setdiff(x2, x1)]

function get_optimal_mapping(rs::RootSystem, vrs::ValRootSystem, metric::Function; optimizer)
    rs, vrs = deepcopy(rs), deepcopy(vrs)
    add_dummy_roots!(rs, vrs)
    return get_optimal_mapping(laterals(rs), laterals(vrs), metric; optimizer)
end

function get_optimal_mapping(rs::Vector{<:Root}, vrs::Vector{<:ValRoot}, metric::Function; optimizer)
    n = length(rs)

    loss_mat = [
        get_loss(rs, vrs, i, j, metric)
        for i in 1:n, j in 1:n
    ]
  
    return get_optimal_mapping(loss_mat; optimizer)
end

function get_optimal_mapping(loss_mat::Matrix{<:Number}; optimizer)
    n = size(loss_mat, 1)
    model = Model(optimizer)
    @variable(model, mm[1:n, 1:n], Bin)
    @objective(model, Min, sum(loss_mat .* mm))
    for i in 1:n
        @constraint(model, sum(mm[i,j] for j in 1:n) == 1)
    end
    for j in 1:n
        @constraint(model, sum(mm[i,j] for i in 1:n) == 1)
    end

    set_silent(model)
    optimize!(model)
    is_solved_and_feasible(model)
    mm = round.(Bool, value.(model[:mm]))
    return [findfirst(c) for c in eachrow(mm)]
end

function get_optimal_mappings(rss::Vector{<:RootSystem}, vrss::Vector{<:ValRootSystem}, metric::Function; optimizer)
    n_og = length(rss)
    rss, vrss = deepcopy(rss), deepcopy(vrss)
    add_dummy_rootsystems!(rss, vrss)
    n = length(rss)

    mapping_mat = [
        get_optimal_mapping(rss[i], vrss[j], metric; optimizer)
        for i in 1:n, j in 1:n
    ]
    loss_mat = [
        get_loss(rss[i], vrss[j], mapping_mat[i, j], metric)
        for i in 1:n, j in 1:n
    ]

    rs_mapping = get_optimal_mapping(loss_mat; optimizer)
    return rs_mapping[1:n_og], [mapping_mat[i, j] for (i, j) in zip(1:n_og, rs_mapping[1:n_og])]
end

function add_dummy_rootsystems!(rss::Vector{<:RootSystem}, vrss::Vector{<:ValRootSystem})
    n_diff = length(rss) - length(vrss)
    n_diff < 0 && add_dummy_rootsystems!(rss, abs(n_diff))
    n_diff > 0 && add_dummy_rootsystems!(vrss, n_diff)

    return nothing
end

function add_dummy_rootsystems!(rss::Vector{RootSystem{T, U}}, n_diff::Integer) where {T, U}
    append!(
        rss, 
        fill(
            RootSystem(
                SimpleRoot(true, RootFragment(false, RootArc{T, U}[])),
                Union{SimpleRoot{T, U}, CompositeRoot{T, U}}[], 
            ),
            n_diff
        )
    )

    return nothing
end

function add_dummy_rootsystems!(vrss::Vector{ValRootSystem{T, U}}, n_diff::Integer) where {T, U}
    append!(
        vrss,
        fill(
            ValRootSystem(
                ValRoot(true, RootEdge{T, U}[], U[], U[]),
                ValRoot{T, U}[]
            ),
            n_diff
        )
    )

    return nothing
end

standard_segments(r::Root) = filter(!is_augmented, og_segment_id.(edges(r)))
standard_segments(vr::ValRoot) = og_segment_id.(edges(vr))

pixels(r::Root) = [(x, y) for (x, y) in zip(xs(r), ys(r))]
pixels(vr::ValRoot) = [(x, y) for (x, y) in zip(xs(vr), ys(vr))]

function plot_mistakes(rg::RootGraph, rss::Vector{<:RootSystem}, vrss::Vector{<:ValRootSystem}, rs_mapping, root_mappings)
    rss, vrss = (deepcopy(rss), deepcopy(vrss))
    for (i, (rs, vrs)) in enumerate(zip(rss, vrss[rs_mapping]))
        add_dummy_roots!(rs, vrs)
        laterals(vrs) .= laterals(vrs)[root_mappings[i]];
    end

    incorrect_segments_ids = [
        [
            full_setdiff(standard_segments(laterals(rs)[i]), standard_segments(laterals(vrs)[i]))
            for i in eachindex(laterals(rs))
        ] |> x -> reduce(vcat, x)
        for (rs, vrs) in zip(rss, vrss[rs_mapping])
    ] |> x -> reduce(vcat, x)

    f = graphplot(rg)
    for re in filter(re -> re.og_segment_id ∈ incorrect_segments_ids, E₀(rg))
        lines!(f[1, 1], rg, re, color = :red, linewidth = 10, alpha = 0.3)
    end
    return f
end