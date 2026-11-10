include("main.jl")
import RootUntangling: read_data, RootEdge, is_primary, edges, primary,
    laterals, roots, V, x, y, xs, ys, rootplot!, rootplot, getsegment, PreGraph,
    Root, SimpleRoot, RootFragment, RootArc, CompositeRoot, RootSystem
import .Makie: HSV

# source

function fix_groundtruth(file::String)
    lines = readlines(file)
    filter!(!isempty, lines)
    split.(lines, ',')
    newtext = [
        length(split(line, ',')) == 2 ? line * ',' : line
        for line in lines
    ] |> x -> *((x .* '\n')...)[1:end-1]

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

    rsa_ids = unique(first.(first.(root_mapping))) |> sort
    root_segment_mappings = [
        [
            mapping[1][2] => mapping[2]
            for mapping in root_mapping
            if mapping[1][1] == rsa_id
        ] |> Dict
        for rsa_id in rsa_ids
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

struct ValRootSystem{T, U}
    primary::ValRoot{T, U}
    laterals::Vector{ValRoot{T, U}}
end
primary(vrs::ValRootSystem) = vrs.primary
laterals(vrs::ValRootSystem) = vrs.laterals
roots(vrs::ValRootSystem) = [primary(vrs); laterals(vrs)]
Base.length(vrs::ValRootSystem) = length(roots(vrs))

function get_rsa_val(rg::RootGraph, segments, root_segment_mapping::Dict)
    primary = ValRoot(
        true,
        filter(e -> e.segment_id ∈ root_segment_mapping[1], E₀(rg)),
        reduce(vcat, xs.(filter(s -> s.id ∈ root_segment_mapping[1], segments))),
        reduce(vcat, ys.(filter(s -> s.id ∈ root_segment_mapping[1], segments))),
    )
    laterals = [
        ValRoot(
            false,
            filter(e -> e.segment_id ∈ root_segment_mapping[i], E₀(rg)),
            reduce(vcat, xs.(filter(s -> s.id ∈ root_segment_mapping[i], segments))),
            reduce(vcat, ys.(filter(s -> s.id ∈ root_segment_mapping[i], segments))),
        )
        for i in filter(x -> x != 1, keys(root_segment_mapping))
    ]
    return ValRootSystem(primary, laterals)
end

function rootplot!(ax::Makie.Axis, vrs::ValRootSystem; kwargs...)
    r = primary(vrs)
    scatter!(ax, r; color = HSV(0, 1, 0), label = "0", markersize = 2, kwargs...)

    for (i, r) in enumerate(laterals(vrs))
        color = HSV(range(0, 360, length = length(vrs))[i], 1, 0.75)
        scatter!(ax, r; color, label = "$i", markersize = 2, kwargs...)
    end
end

function rootplot(vrs::ValRootSystem; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, vrs; line_kwargs...)

    return f
end

function add_dummy_roots!(rsa::RootSystem, rsa_val::ValRootSystem)
    laterals_pred = laterals(rsa)
    laterals_val = laterals(rsa_val)
    n_diff = length(laterals_pred) - length(laterals_val)
    n_diff < 0 && add_dummy_roots!(laterals_pred, abs(n_diff))
    n_diff > 0 && add_dummy_roots!(laterals_val, n_diff)

    return nothing
end

function add_dummy_roots!(roots::Vector{Union{CompositeRoot{T, U}, SimpleRoot{T, U}}}, n_diff::Integer) where {T, U}
    append!(roots, fill(SimpleRoot(false, RootFragment(false, RootArc{T}[])), n_diff))

    return nothing
end

function add_dummy_roots!(roots::Vector{ValRoot{T, U}}, n_diff::Integer) where {T, U}
    append!(roots, fill(ValRoot(false, RootEdge{T, U}[], U[], U[]), n_diff))

    return nothing
end

function get_loss(rsa::RootSystem, rsa_val::ValRootSystem, root_mapping::Vector{<:Integer}, metric::Function;
        verbose::Bool = false
    )
    # primary root
    primary_dist = get_loss(primary(rsa), primary(rsa_val), metric)

    # laterals
    lateral_dist = 0
    for i in eachindex(root_mapping)
        lateral_dist += get_loss(laterals(rsa), laterals(rsa_val), i, root_mapping[i], metric)
        verbose && @info "Loss in root $(i): $(get_loss(laterals(rsa), laterals(rsa_val), i, root_mapping[i], metric))"
    end

    return primary_dist + lateral_dist
end

get_loss(rs::Vector{<:Root}, rvs::Vector{<:ValRoot}, i::Integer, j::Integer, metric::Function) = get_loss(rs[i], rvs[j], metric)

function get_loss(r::Root, rv::ValRoot, metric::Function)
    return length(full_setdiff(metric(r), metric(rv)))
end

full_setdiff(x1::AbstractArray, x2::AbstractArray) = [setdiff(x1, x2); setdiff(x2, x1)]

function get_optimal_mapping(rs::Vector{<:Root}, rsv::Vector{<:ValRoot}, metric::Function; optimizer)
    @assert length(rs) == length(rsv)
    n = length(rs)

    loss_mat = [
        get_loss(rs, rsv, i, j, metric)
        for i in 1:n, j in 1:n
    ]
  
    model = Model(optimizer)
    @variable(model, mm[1:n, 1:n], Bin)
    @objective(model, Min, sum(loss_mat .* mm))
    for i in 1:n
        @constraint(model, sum(mm[i,j] for j in 1:n) == 1)
    end
    for j in 1:n
        @constraint(model, sum(mm[i,j] for i in 1:n) == 1)
    end

    optimize!(model)
    is_solved_and_feasible(model)
    mm = round.(Bool, value.(model[:mm]))
    return [findfirst(r) for r in eachrow(mm)]
end

standard_segments(r::Root) = filter(!is_augmented, segment_id.(edges(r)))
standard_segments(vr::ValRoot) = segment_id.(edges(vr))

pixels(r::Root) = [(x, y) for (x, y) in zip(xs(r), ys(r))]
pixels(vr::ValRoot) = [(x, y) for (x, y) in zip(xs(vr), ys(vr))]

# script

# categories = ["baby", "easy", "intermediate", "tough"]
# n_rois = [3, 6, 6, 4]
# for (i, category) in enumerate(categories), roi in 1:n_rois[i]
#     file = "data/validation/$(category)/ROI_$(roi)/GT.txt"
#     fix_groundtruth(file)
# end

begin
    file = "data/validation/$(difficulty)/ROI_$(roi_nr)/GT.txt"

    valdata = RootUntangling.read_data(file, :Segment_ID, missing_string = "-", parse_strings = false) |>
    x -> [parse(Int64, key) => x[key] for key in keys(x)] |> Dict
    root_segment_mappings = get_root_segment_mappings(valdata);
    root_segmentation_errors = [
        filter(!ismissing, Any[
            valdata[seg_id][:real_length]
            for seg_id in reduce(vcat, values(root_segment_mapping))
        ]) .|> (x -> parse(Int64, x)) |> sum
        for root_segment_mapping in root_segment_mappings
    ]
end

rsa = rss_new[1]

begin
    vertex_data_dict, edge_data_dict = RootUntangling.get_data_dicts(filename_segments, filename_vertices; reverse_y, angle_correction)
    metavertices = RootUntangling.getmetavertices(vertex_data_dict)
    segments = RootUntangling.getsegments(metavertices, vertex_data_dict, edge_data_dict)
    RootUntangling.sort_segment_coordinates!(segments, metavertices)

    rsa_nr = 3
    rsa_val = get_rsa_val(rg_full, segments, root_segment_mappings[rsa_nr])
    f = Figure()
    rootplot!(Axis(f[1, 1]), rsa)
    f[1, 2] = Legend(f, f.content[end])
    rootplot!(Axis(f[1, 3]), rsa_val)
    f[1, 4] = Legend(f, f.content[end])
    f
end

add_dummy_roots!(rsa, rsa_val)

using Random
n = length(laterals(rsa))
metric = pixels

losses = [get_loss(rsa, rsa_val, shuffle(1:n), metric) for _ in 1:1_000]
hist(losses)
minimum(losses)

optimal_mapping = get_optimal_mapping(laterals(rsa), laterals(rsa_val), metric, optimizer = HiGHS.Optimizer)
laterals(rsa_val) .= laterals(rsa_val)[optimal_mapping]
final_loss = get_loss(rsa, rsa_val, collect(1:n), metric, verbose = true)

begin
    f = Figure()
    rootplot!(Axis(f[1, 1]), rsa_val)
    f[1, 2] = Legend(f, f.content[end])
    rootplot!(Axis(f[1, 3]), rsa)
    f[1, 4] = Legend(f, f.content[end])
    Axis(f[1, 1]).title = "GT"
    Axis(f[1, 3]).title = "Predicted (loss = $(final_loss) pixels)"
    f
end

# investigations

baddies = reduce(vcat, [getsegments(pg, i) for i in [263, 264, 508, 509]])

f = Figure(size = (600, 600))
lines(f[1, 1], xs(baddies[1]), ys(baddies[1]))
for i in 2:length(baddies)
    lines!(f[1, 1], xs(baddies[i]), ys(baddies[i]))
end
f