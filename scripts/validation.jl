include("main.jl")
import RootUntangling: read_data, RootEdge, is_primary, edges, primary,
    laterals, roots, V, x, y, xs, ys, rootplot!, rootplot,
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
end
is_primary(vr::ValRoot) = vr.is_primary
edges(vr::ValRoot) = vr.edges
Base.length(vr::ValRoot) = length(edges(vr))

xs(rg::RootGraph, vr::ValRoot) = reduce(vcat, [
    [x.(V(rg, e)); NaN]
    for e in vr.edges
], init = Float64[])
ys(rg::RootGraph, vr::ValRoot) = reduce(vcat, [
    [y.(V(rg, e)); NaN]
    for e in vr.edges
], init = Float64[])
Makie.convert_arguments(::Type{<:Lines}, rg::RootGraph, vr::ValRoot) = (xs(rg, vr), ys(rg, vr))

struct ValRootSystem{T, U}
    primary::ValRoot{T, U}
    laterals::Vector{ValRoot{T, U}}
end
primary(vrs::ValRootSystem) = vrs.primary
laterals(vrs::ValRootSystem) = vrs.laterals
roots(vrs::ValRootSystem) = [primary(vrs); laterals(vrs)]
Base.length(vrs::ValRootSystem) = length(roots(vrs))

function get_rsa_val(rg::RootGraph, root_segment_mapping::Dict)
    primary = ValRoot(true, filter(e -> e.segment_id ∈ root_segment_mapping[1], E₀(rg)))
    laterals = [
        ValRoot(false, filter(e -> e.segment_id ∈ root_segment_mapping[i], E₀(rg)))
        for i in filter(x -> x != 1, keys(root_segment_mapping))
    ]
    return ValRootSystem(primary, laterals)
end

function rootplot!(ax::Makie.Axis, rg::RootGraph, vrs::ValRootSystem; kwargs...)
    r = primary(vrs)
    lines!(ax, rg, r; color = HSV(0, 1, 0), linestyle = :solid, label = "0", linewidth = 2, kwargs...)

    for (i, r) in enumerate(laterals(vrs))
        color = HSV(range(0, 360, length = length(vrs))[i], 1, 0.75)
        lines!(ax, rg, r; color, linestyle = :dot, label = "$i", linewidth = 2, kwargs...)
    end
end

function rootplot(rg::RootGraph, vrs::ValRootSystem; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, rg, vrs; line_kwargs...)

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
    append!(roots, fill(ValRoot(false, RootEdge{T, U}[]), n_diff))

    return nothing
end

function get_loss(rsa::RootSystem, rsa_val::ValRootSystem, root_mapping::Vector{<:Integer}; 
        verbose::Bool = false
    )
    # primary root
    primary_dist = get_loss(primary(rsa), primary(rsa_val))

    # laterals
    lateral_dist = 0
    for i in eachindex(root_mapping)
        lateral_dist += get_loss(laterals(rsa), laterals(rsa_val), i, root_mapping[i])
        verbose && @info "Loss in root $(i): $(get_loss(laterals(rsa), laterals(rsa_val), i, root_mapping[i]))"
    end

    return primary_dist + lateral_dist
end

get_loss(rs::Vector{<:Root}, rvs::Vector{<:ValRoot}, i::Integer, j::Integer) = get_loss(rs[i], rvs[j])

function get_loss(r::Root, rv::ValRoot)
    segments = r |> edges .|> segment_id
    filter!(!is_augmented, segments)
    segments_val = rv |> edges .|> segment_id
    return length(full_setdiff(segments, segments_val)) #! needs to become pixels
end

full_setdiff(x1::AbstractArray, x2::AbstractArray) = [setdiff(x1, x2); setdiff(x2, x1)]

function get_optimal_mapping(rs::Vector{<:Root}, rsv::Vector{<:ValRoot}; optimizer)
    @assert length(rs) == length(rsv)
    n = length(rs)

    loss_mat = [
        get_loss(rs, rsv, i, j)
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
    rsa_nr = 2
    rsa_val = get_rsa_val(rg_full, root_segment_mappings[rsa_nr])
    f = Figure()
    rootplot!(Axis(f[1, 1]), rg, rsa)
    f[1, 2] = Legend(f, f.content[end])
    rootplot!(Axis(f[1, 3]), rg_full, rsa_val)
    f[1, 4] = Legend(f, f.content[end])
    f
end

add_dummy_roots!(rsa, rsa_val)

using Random
n = length(laterals(rsa))
losses = [get_loss(rsa, rsa_val, shuffle(1:n)) for _ in 1:100_000]
hist(losses, bins = 0:25)
minimum(losses)

optimal_mapping = get_optimal_mapping(laterals(rsa), laterals(rsa_val), optimizer = HiGHS.Optimizer)
laterals(rsa_val) .= laterals(rsa_val)[optimal_mapping]
get_loss(rsa, rsa_val, collect(1:n), verbose = true)

begin
    f = Figure()
    rootplot!(Axis(f[1, 1]), rg, rsa)
    f[1, 2] = Legend(f, f.content[end])
    rootplot!(Axis(f[1, 3]), rg_full, rsa_val)
    f[1, 4] = Legend(f, f.content[end])
    f
end

# investigations

bastard_id = 227

f = graphplot(rg_full, augmented_alpha = 0)
bastard = E₀(rg_full)[findfirst(e -> e.segment_id == bastard_id, E₀(rg_full))]
lines!(f.content[1], rg_full, bastard, color = "red", linewidth = 10)
f