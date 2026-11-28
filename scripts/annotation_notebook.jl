### A Pluto.jl notebook ###
# v1.0.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ bace35f3-a1e4-49ff-9f63-fb4ca65d2176
using Pkg

# ╔═╡ 7d5639a6-a788-4877-8ac8-394ffb45a220
begin
	Pkg.activate()
	using Revise, PlutoUI
end

# ╔═╡ 1ebe1d3c-dd67-44d5-b292-c6ce4b9f9aa1
begin
	Pkg.activate(".")
	using RootUntangling
	using WGLMakie
	using JuMP, HiGHS
	import .Makie.HSV, .Makie.HSVA, RootUntangling.mean
end

# ╔═╡ 1423c220-fb6d-4c4c-bd07-140e29cc79cc
md"## Imports"

# ╔═╡ 71e10643-4c1d-4069-bd9e-eb8fd58dccc7
md"## Data reading"

# ╔═╡ 72b6f5ea-f063-4705-bc9c-b06d55fd0d43
begin
    cluster_threshold = 3
    loop_threshold = 5
    reverse_y = true
    angle_correction = -
end;

# ╔═╡ 9bedbeba-0a6c-447c-ae17-0cec8e95d23c
@bind dirname Select(["playground", "validation", "calibration"])

# ╔═╡ 079900fc-3ce2-4b67-9de1-3f2654523eff
if dirname == "validation"
	@bind difficulty Select(["baby", "easy", "intermediate", "tough", "ungodly"])
end

# ╔═╡ 8d59a0ee-4f91-48b7-993c-20e799c957ce
begin
    min_vertices = 5
    y_min = -500
    y_max = 5000
end

# ╔═╡ 00312cba-17c3-4bf7-acf2-e637a99c710e
md"## Annotation"

# ╔═╡ d16d027b-f57b-4d91-bafd-4ea0d14241f6
@bindname from_scratch Switch()

# ╔═╡ fee7a9e4-481e-40c9-9ffc-fa12fcf1b12b
if !from_scratch
	@bindname load_corrected Switch()
end

# ╔═╡ af73fbe5-f5ba-494e-a82a-7391eab778c0
if from_scratch
	@bind go CounterButton("Run model")
end

# ╔═╡ c9f608c8-0088-4667-909e-2ed0cc0d1866
@bind done CounterButton("Save annotation")

# ╔═╡ 2e87558f-1c4f-40d7-b98f-8aa4e8858bb8
md"## Just checking some stuff aha"

# ╔═╡ 150e2c5e-9828-4180-98b1-cb027d76cba6
md"## Source"

# ╔═╡ 92477ec2-047e-4d6b-971a-9985662de553
function get_ROIs(dir::String)
	names = readdir("../data/$(dir)")
	filter!(name -> occursin("ROI", name), names)
	return parse.(Int64, last.(split.(names, '_')))
end

# ╔═╡ 6c3bd7e0-1d81-46e0-935d-4f6ec23a55a2
function input_cluster_root_counts(rgs::Vector{<:RootGraph})
	return PlutoUI.combine() do Child
		
		inputs = [
			md"""Cluster $(i): $(
				Child(NumberField(default = 1))
			)"""
			
			for i in eachindex(rgs)
		]
		
		md"""
		#### Root counts per cluster
		$(inputs)
		"""
	end
end

# ╔═╡ ff4f7aee-13c6-49c2-837d-12a3e741c28c
root_color(root_idx, num_roots) = (
    root_idx == 1 ? HSV(0, 1, 0) : 
    num_roots == 1 ? HSV(0, 1, 0.75) : 
    root_idx <= num_roots ? HSV(range(0, 360, length = num_roots)[root_idx], 1, 0.75) : HSV(root_idx * 360 * (1 + sqrt(5))/2, 1, 1)
) #!

# ╔═╡ 16b156b2-a73b-424c-b271-9e15994c4d2f
function set_colors!(edge_color_dict::Dict, annotation_dict::Dict, 
		rg::RootGraph, root_nums::Vector
	)
	for re in E₀(rg)
		haskey(annotation_dict, re.og_segment_id) || continue
		colors = [
			root_color(root_idx, root_nums[rs_idx])
			for (rs_idx, root_idx) in annotation_dict[re.og_segment_id]
		]
		if !isempty(colors)
			edge_color_dict[re][] = Makie.weighted_color_mean(ones(length(colors))/length(colors), colors)
		else
			edge_color_dict[re][] = HSV(0, 0, 0.5)
		end
	end

	return nothing
end

# ╔═╡ 1de98531-6b6a-4443-8fb7-793509730fe8
function set_annotations!(annotations::Observable, annotation_dict::Dict, 
		rg::RootGraph
	)
	
	annotations[] = [
		string(
			get(annotation_dict, re.og_segment_id, "[]") |>
				x -> isempty(x) ? "[]" : x
		)
		for re in E₀(rg)
	]

	return nothing
end

# ╔═╡ 7f2bd2eb-83d1-495e-b4bf-226beca9f733
function get_segment_id(f::Makie.Figure)
	elements = Makie.pick_sorted(f.scene, events(f).mouseposition[], 30)
    line_idx = findfirst(x -> x[1] isa Lines, elements)
    isnothing(line_idx) && return nothing
    
    line_element = elements[line_idx]
    re_picked = line_element[1].arg2.value[]
    return re_picked.og_segment_id
end

# ╔═╡ 75f65731-fcd7-4d60-baee-3b26f9a3a95d
if @isdefined difficulty
	dir = dirname * '/' * difficulty
else
	dir = dirname
end

# ╔═╡ 690ea173-af44-475e-bc5a-1f1044e48309
@bind roi_nr Select(get_ROIs(dir))

# ╔═╡ d88b834b-241f-4972-ae13-281d8b62b26d
begin
	filename_segments = "../data/$(dir)/ROI_$(roi_nr)/segment_info_with_coords.csv"
    filename_vertices = "../data/$(dir)/ROI_$(roi_nr)/bp1_segments_grouped.csv"
	filename_GT = "../data/$(dir)/ROI_$(roi_nr)/$(@isdefined(load_corrected) && load_corrected ? "GT_new.txt" : "GT.txt")"
end

# ╔═╡ c26c9948-8ff1-44c3-8b68-8e0907fdb3ec
begin
    pg = get_pregraph(filename_segments, filename_vertices; cluster_threshold, loop_threshold, reverse_y, angle_correction)
    rg_full = get_rootgraph(pg)
end

# ╔═╡ 792c91d6-4863-4165-b926-c8b0769a0481
begin
    rgs = get_subgraphs(rg_full) |>
        rgs -> filter(rg -> length(rg) > min_vertices, rgs) |>
        rgs -> filter(rg -> minimum(y.(V₀(rg))) >= y_min && maximum(y.(V₀(rg))) <= y_max, rgs) |>
        rgs -> sort(rgs, by = rg -> mean(x.(V₀(rg))));
    
    f_multi = Figure()
    ax_multi = Axis(f_multi[1, 1]; aspect = DataAspect())
    
    ls = [
        lines!(ax_multi, rg, E₀(rg), color = Makie.HSV(i / length(rgs) * 360, 1, 1))
        for (i, rg) in enumerate(rgs)
    ]
    Legend(f_multi[1, 2], ls, string.(1:length(ls)))

    f_multi
end

# ╔═╡ 2676da59-c62a-4595-9d05-7a36389a0196
if from_scratch
	@bind cluster_root_counts_input input_cluster_root_counts(rgs)
end

# ╔═╡ ca219486-3522-4f3a-837b-b3fc99b2e670
if !from_scratch
    done
    vrss = get_valrootsystems(filename_GT, pg, rg_full)
    rootplot(vrss, size = (1000, 500))
end |> WideCell

# ╔═╡ 0a1e3cc5-d1bb-4049-b4a3-79b05e46bcd4
if from_scratch
	if go > 0
		cluster_root_counts = [cluster_root_counts_input...]
		root_clusters = RootUntangling.get_root_clusters(
			rgs, pg, cluster_root_counts; optimizer = HiGHS.Optimizer, time_limit = 60
		)
		annotation_dict = RootUntangling.get_annotation_dict(rgs, root_clusters)
		error_dict = Dict()
	end
else
	if isfile(filename_GT)
		RootUntangling.fix_groundtruth(filename_GT)
		annotation_dict = read_annotation_dict(filename_GT)
		error_dict = read_error_dict(filename_GT)
	end
end

# ╔═╡ 17648059-dcb1-400a-bc89-1a5bb734ee03
if @isdefined annotation_dict
    root_nums = [
        length(unique(last.(filter(k -> k[1] == i, reduce(vcat, values(annotation_dict))))))
        for i in 1:maximum(first.(reduce(vcat, values(annotation_dict))))
    ]
    # root_nums = fill(100, 11)
else
    @warn "There is no GT file available for this root system"
end

# ╔═╡ 08de3104-630e-4792-836f-0ab11d324470
if @isdefined annotation_dict
	@bind working_on_button Select(
		unique(reduce(vcat, values(annotation_dict))) |>
			x -> [x; [(rs, maximum(last.(x[findall(_x -> _x[1] == rs, x)]))+1) for rs in unique(first.(x))]] |>
			x -> [x; [(maximum(unique(first.(x)))+1, 1)]] |>
			sort
	)
end

# ╔═╡ 2a4561ca-f05e-4fba-b14a-2d5ea5de9f96
if @isdefined annotation_dict
    # instantiate figure with observables
    f = Figure(; size = (1600, 900))
    ax = Axis(f[1, 1]; aspect = DataAspect())

    edge_color_dict = [
        re => Observable(HSVA(0, 0, 0, 0.3))
        for re in E₀(rg_full)
    ] |> Dict

    annotations = [
        ""
        for re in E₀(rg_full)
    ] |> Observable

    for re in E₀(rg_full)
        lines!(ax, rg_full, re; color = edge_color_dict[re], linestyle = :solid)
    end

    x_coords = [mean(xs(rg_full, re)) for re in E₀(rg_full)]
    y_coords = [mean(ys(rg_full, re)) for re in E₀(rg_full)]
    text!(ax, x_coords, y_coords; text = annotations, fontsize = 8)

    set_colors!(edge_color_dict, annotation_dict, rg_full, root_nums)
    set_annotations!(annotations, annotation_dict, rg_full)

    info_label = Observable("(1, 1)")
    text!(
        ax, 0.05, 0.95;
        text = info_label,
        space = :relative,
        fontsize = 16
    )
    debug_label = Observable("0")
    text!(
        ax, 0.95, 0.95;
        text = debug_label,
        space = :relative,
        fontsize = 16
    )
    
    # update based on user input
    s_id = nothing
    working_on = Observable((1, 1))

    on(events(ax.scene).mouseposition) do event
        s_id = get_segment_id(f)
        debug_label[] = string(s_id)
    end
    
    on(events(f).mousebutton) do event
        if event.button == Mouse.left && event.action == Mouse.press
            s_id = get_segment_id(f)
            if !isnothing(s_id)
                current_annotation_idx = findfirst(
                    x -> x == working_on[], 
                    get!(annotation_dict, s_id, Tuple{Int64, Int64}[])
                )
    
                if isnothing(current_annotation_idx)
                    push!(annotation_dict[s_id], working_on[])
                else
                    deleteat!(annotation_dict[s_id], current_annotation_idx)
                end
                set_colors!(edge_color_dict, annotation_dict, rg_full, root_nums)
                set_annotations!(annotations, annotation_dict, rg_full)
            end
        end
    end

    on(events(f).keyboardbutton) do event
        if event.action == Keyboard.press
            if event.key == Keyboard.space
                s_id = get_segment_id(f)
                if haskey(annotation_dict, s_id) && !isempty(annotation_dict[s_id])
                    w_idx = findfirst(id -> id == working_on[], annotation_dict[s_id])
                    working_on[] = isnothing(w_idx) ? annotation_dict[s_id][1] : 
                        get(annotation_dict[s_id], w_idx + 1, annotation_dict[s_id][1])
                    info_label[] = string(working_on[])
                end
            elseif event.key == Keyboard.q
                working_on[] = (working_on[][1], working_on[][2] - 1)
                info_label[] = string(working_on[])
            elseif event.key == Keyboard.w
                working_on[] = (working_on[][1], working_on[][2] + 1)
                info_label[] = string(working_on[])
                info_label[] = string(working_on[])
            elseif event.key == Keyboard.e
                working_on[] = reduce(vcat, values(annotation_dict)) |>
                    unique |>
                    x -> filter(i -> i[1] == working_on[][1], x) |>
                    x -> (working_on[][1], maximum(last.(x)) + 1)
                info_label[] = string(working_on[])
            end
        end
    end
    f
end |> x -> WideCell(x; max_width = 1600)

# ╔═╡ 94418555-83fb-4610-a04d-89d2823d3efb
[k => v for (k, v) in error_dict if !ismissing(v)]

# ╔═╡ 1053aa28-910d-49b0-a85b-ab01d2610f9c
if @isdefined root_clusters
	metric = pixels
    rss = reduce(vcat, root_clusters)
	rs_mapping, root_mappings = get_optimal_mappings(rss, vrss, metric, optimizer = HiGHS.Optimizer)

	for (i, rs) in enumerate(rss)
        vrs = vrss[rs_mapping[i]]
        rs_copy, vrs_copy = deepcopy(rs), deepcopy(vrs);
        add_dummy_roots!(rs_copy, vrs_copy)
        final_loss = get_loss(rs_copy, vrs_copy, root_mappings[i], metric, verbose = true)
        laterals(vrs_copy) .= laterals(vrs_copy)[root_mappings[i]];
    
        local f = rootplot(vrs_copy, rs_copy)
        Axis(f[1, 1]).title = "GT"
        Axis(f[1, 3]).title = "Predicted"
        Label(f[0, :], "loss = $(final_loss) $(string(metric))")
        WGLMakie.Bonito.export_static(
            homedir() * "/Downloads/validation_$(difficulty)_$(roi_nr)_$i.html",
            WGLMakie.Bonito.App(f)
        )
    end
end

# ╔═╡ e756a41b-cf15-4b3a-bfe7-3cb4ddde0ac4
if @isdefined working_on
	working_on[] = working_on_button
end

# ╔═╡ 8f725f54-4429-48a4-856f-93d9b664b74e
done_ref = Ref(done)

# ╔═╡ 78a5a17f-8064-437f-bdf5-6ac47fe0741c
if done_ref[] > 0
    write_annotation("../data/$(dir)/ROI_$(roi_nr)/GT_new.txt", rg_full, annotation_dict, error_dict)
    done_ref[] = 0
    println("Annotation saved$("!"^(ceil(Int64, 30*rand())))")
end

# ╔═╡ Cell order:
# ╟─1423c220-fb6d-4c4c-bd07-140e29cc79cc
# ╠═bace35f3-a1e4-49ff-9f63-fb4ca65d2176
# ╠═7d5639a6-a788-4877-8ac8-394ffb45a220
# ╠═1ebe1d3c-dd67-44d5-b292-c6ce4b9f9aa1
# ╟─71e10643-4c1d-4069-bd9e-eb8fd58dccc7
# ╠═72b6f5ea-f063-4705-bc9c-b06d55fd0d43
# ╟─9bedbeba-0a6c-447c-ae17-0cec8e95d23c
# ╟─079900fc-3ce2-4b67-9de1-3f2654523eff
# ╟─690ea173-af44-475e-bc5a-1f1044e48309
# ╠═8d59a0ee-4f91-48b7-993c-20e799c957ce
# ╟─792c91d6-4863-4165-b926-c8b0769a0481
# ╟─00312cba-17c3-4bf7-acf2-e637a99c710e
# ╟─17648059-dcb1-400a-bc89-1a5bb734ee03
# ╟─d16d027b-f57b-4d91-bafd-4ea0d14241f6
# ╟─fee7a9e4-481e-40c9-9ffc-fa12fcf1b12b
# ╟─2676da59-c62a-4595-9d05-7a36389a0196
# ╟─af73fbe5-f5ba-494e-a82a-7391eab778c0
# ╟─08de3104-630e-4792-836f-0ab11d324470
# ╟─2a4561ca-f05e-4fba-b14a-2d5ea5de9f96
# ╟─c9f608c8-0088-4667-909e-2ed0cc0d1866
# ╟─78a5a17f-8064-437f-bdf5-6ac47fe0741c
# ╟─2e87558f-1c4f-40d7-b98f-8aa4e8858bb8
# ╠═ca219486-3522-4f3a-837b-b3fc99b2e670
# ╠═94418555-83fb-4610-a04d-89d2823d3efb
# ╟─1053aa28-910d-49b0-a85b-ab01d2610f9c
# ╟─150e2c5e-9828-4180-98b1-cb027d76cba6
# ╟─92477ec2-047e-4d6b-971a-9985662de553
# ╟─6c3bd7e0-1d81-46e0-935d-4f6ec23a55a2
# ╟─ff4f7aee-13c6-49c2-837d-12a3e741c28c
# ╟─16b156b2-a73b-424c-b271-9e15994c4d2f
# ╟─1de98531-6b6a-4443-8fb7-793509730fe8
# ╟─7f2bd2eb-83d1-495e-b4bf-226beca9f733
# ╟─d88b834b-241f-4972-ae13-281d8b62b26d
# ╟─c26c9948-8ff1-44c3-8b68-8e0907fdb3ec
# ╟─75f65731-fcd7-4d60-baee-3b26f9a3a95d
# ╟─0a1e3cc5-d1bb-4049-b4a3-79b05e46bcd4
# ╟─e756a41b-cf15-4b3a-bfe7-3cb4ddde0ac4
# ╟─8f725f54-4429-48a4-856f-93d9b664b74e
