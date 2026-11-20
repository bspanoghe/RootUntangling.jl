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
	using Revise, PlutoUI; TableOfContents()
end

# ╔═╡ 1ebe1d3c-dd67-44d5-b292-c6ce4b9f9aa1
begin
	Pkg.activate(".")
	using RootUntangling
	using WGLMakie
	using JuMP, HiGHS
	import .Makie.HSV, RootUntangling.mean
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

# ╔═╡ 0bce83d3-e008-45d0-b32c-0b08c9b6711b
md"### Rootsystem selection"

# ╔═╡ 8d59a0ee-4f91-48b7-993c-20e799c957ce
begin
    min_vertices = 20
    y_min = -1000
    y_max = 4800
end

# ╔═╡ 3fc9f8b8-e6e8-4409-a920-98d30ac9c2d1
md"## Data visualisation"

# ╔═╡ 690ea173-af44-475e-bc5a-1f1044e48309
@bind roi_nr Select(1:7, default = 4)

# ╔═╡ c26c9948-8ff1-44c3-8b68-8e0907fdb3ec
begin
    filename_segments = "../data/playground/ROI_$(roi_nr)/segment_info_with_coords.csv"
    filename_vertices = "../data/playground/ROI_$(roi_nr)/bp1_segments_grouped.csv"
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

# ╔═╡ d86ed7ae-877e-42e8-a748-109e92daab38
@bind subidx Select(1:length(rgs), default = 1)

# ╔═╡ 3f40ee53-68ad-4bfa-a784-bbad4a366f02
begin
	rg = rgs[subidx]
	graphplot(rg, augmented_alpha = 0.01)
end

# ╔═╡ 7894d17f-b98c-4e19-8d3b-4fb690aaeb42
md"## Solving"

# ╔═╡ 990ad7fb-4256-41ce-a6fe-520a103c5ea7
begin
	time_limit = 60
	num_roots = 1
	max_overlapping = 5
end;

# ╔═╡ 5355a78b-3bee-454e-a18b-0e6ac2d202b7
# ╠═╡ show_logs = false
begin
    model, time = @timed solve_rsa(
        rg; optimizer = HiGHS.Optimizer, time_limit, num_roots, max_overlapping
    );
    rss_raw = get_rootsystems(rg, pg, model);
    rss = greedy_switch(rss_raw, f_obj = weighted_tortuosity);
end;

# ╔═╡ 20ffd4b3-481b-4e7b-83d2-7be41bd497f4
begin
	f_root = Figure()
	RootUntangling.rootplot!(
		Axis(f_root[1, 1], aspect = DataAspect()), rss_raw, size = (600, 600), 
		title = "Weighted tortuosity: $(weighted_tortuosity(rss_raw))"
	)
	RootUntangling.rootplot!(
		Axis(f_root[1, 2], aspect = DataAspect()), rss, size = (600, 600),
		title = "Weighted tortuosity: $(weighted_tortuosity(rss_raw))"
	)
	f_root
end

# ╔═╡ 00312cba-17c3-4bf7-acf2-e637a99c710e
md"## Annotation"

# ╔═╡ f7d8603d-9939-4b79-bc1b-27a3fc311234
# ╠═╡ show_logs = false
root_clusters = RootUntangling.get_root_clusters(rgs, pg, [1, 1]; optimizer = HiGHS.Optimizer, time_limit = 60);

# ╔═╡ 95fce8b3-1db5-4a61-b4e3-0e8e36181416
annotation_dict = RootUntangling.get_annotation_dict(rgs, root_clusters);

# ╔═╡ 17648059-dcb1-400a-bc89-1a5bb734ee03
root_nums = [
    length(unique(last.(filter(k -> k[1] == i, reduce(vcat, values(annotation_dict))))))
    for i in sort(unique(first.(reduce(vcat, values(annotation_dict)))))
]

# ╔═╡ ff4f7aee-13c6-49c2-837d-12a3e741c28c
root_color(root_idx, num_roots) = (
    root_idx == 1 ? HSV(0, 1, 0) : (
        num_roots == 1 ? HSV(0, 1, 0.75) : 
            HSV(range(0, 360, length = num_roots)[root_idx], 1, 0.75)
    )
)

# ╔═╡ 16b156b2-a73b-424c-b271-9e15994c4d2f
function set_colors!(edge_color_dict::Dict, annotation_dict::Dict, 
		rg::RootGraph, root_nums::Vector
	)
	for re in E₀(rg)
		colors = [
			root_color(root_idx, root_nums[rs_idx])
			for (rs_idx, root_idx) in annotation_dict[segment_id(re)]
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
function set_annotations!(edge_annotation_dict::Dict, annotation_dict::Dict, 
		rg::RootGraph
	)
	for re in E₀(rg)
		ann = annotation_dict[segment_id(re)]
		edge_annotation_dict[re][] = isempty(ann) ? "[]" : string(ann)
	end

	return nothing
end

# ╔═╡ 08de3104-630e-4792-836f-0ab11d324470
@bind working_on Select(
	sort(unique(reduce(vcat, values(annotation_dict))))
)

# ╔═╡ 2a4561ca-f05e-4fba-b14a-2d5ea5de9f96
begin    
    f = Figure(; size = (500, 1000))
    ax = Axis(f[1, 1]; aspect = DataAspect())

    edge_color_dict = [
        re => Observable(HSV(0, 1, 0))
        for re in E₀(rg)
    ] |> Dict

    edge_annotation_dict = [
        re => Observable("")
        for re in E₀(rg)
    ] |> Dict

    for re in E₀(rg)
        lines!(ax, rg, re; color = edge_color_dict[re], linestyle = :solid)
        annotation!(
            f.content[1], (mean(xs(rg, re)), mean(ys(rg, re))); 
            text = edge_annotation_dict[re], color = :black, 
            shrink = (0, 0), fontsize = 8
        )
    end

    set_colors!(edge_color_dict, annotation_dict, rg, root_nums)
    set_annotations!(edge_annotation_dict, annotation_dict, rg)

    on(events(f).mousebutton, priority = 2) do event
        if event.button == Mouse.left && event.action == Mouse.press
            elements = Makie.pick_sorted(f.scene, events(f).mouseposition[], 30)
            line_idx = findfirst(x -> x[1] isa Lines, elements)
            if !isnothing(line_idx)
                line_element = elements[line_idx]
                re_picked = line_element[1].arg2.value[]
                
                s_id = segment_id(re_picked)

                current_annotation_idx = findfirst(
                    x -> x == working_on, annotation_dict[s_id]
                )

                if isnothing(current_annotation_idx)
                    push!(annotation_dict[s_id], working_on)
                else
                    deleteat!(annotation_dict[s_id], current_annotation_idx)
                end
                set_colors!(edge_color_dict, annotation_dict, rg, root_nums)
                set_annotations!(edge_annotation_dict, annotation_dict, rg)
            end
        end
    end
    f
end

# ╔═╡ 9b60de2c-19e6-40ca-b4b0-bc60b19f3adb
annotation_dict

# ╔═╡ c9f608c8-0088-4667-909e-2ed0cc0d1866
@bind done CounterButton()

# ╔═╡ 78a5a17f-8064-437f-bdf5-6ac47fe0741c
if done > 0
    write_annotation(homedir() * "/Downloads/test_$(roi_nr).txt", rg_full, rgs, root_clusters)
end

# ╔═╡ Cell order:
# ╟─1423c220-fb6d-4c4c-bd07-140e29cc79cc
# ╠═bace35f3-a1e4-49ff-9f63-fb4ca65d2176
# ╠═7d5639a6-a788-4877-8ac8-394ffb45a220
# ╠═1ebe1d3c-dd67-44d5-b292-c6ce4b9f9aa1
# ╟─71e10643-4c1d-4069-bd9e-eb8fd58dccc7
# ╠═72b6f5ea-f063-4705-bc9c-b06d55fd0d43
# ╠═c26c9948-8ff1-44c3-8b68-8e0907fdb3ec
# ╟─0bce83d3-e008-45d0-b32c-0b08c9b6711b
# ╠═8d59a0ee-4f91-48b7-993c-20e799c957ce
# ╠═792c91d6-4863-4165-b926-c8b0769a0481
# ╟─3fc9f8b8-e6e8-4409-a920-98d30ac9c2d1
# ╠═690ea173-af44-475e-bc5a-1f1044e48309
# ╠═d86ed7ae-877e-42e8-a748-109e92daab38
# ╟─3f40ee53-68ad-4bfa-a784-bbad4a366f02
# ╟─7894d17f-b98c-4e19-8d3b-4fb690aaeb42
# ╠═990ad7fb-4256-41ce-a6fe-520a103c5ea7
# ╟─5355a78b-3bee-454e-a18b-0e6ac2d202b7
# ╟─20ffd4b3-481b-4e7b-83d2-7be41bd497f4
# ╟─00312cba-17c3-4bf7-acf2-e637a99c710e
# ╠═f7d8603d-9939-4b79-bc1b-27a3fc311234
# ╠═95fce8b3-1db5-4a61-b4e3-0e8e36181416
# ╟─17648059-dcb1-400a-bc89-1a5bb734ee03
# ╟─ff4f7aee-13c6-49c2-837d-12a3e741c28c
# ╟─16b156b2-a73b-424c-b271-9e15994c4d2f
# ╟─1de98531-6b6a-4443-8fb7-793509730fe8
# ╟─08de3104-630e-4792-836f-0ab11d324470
# ╠═2a4561ca-f05e-4fba-b14a-2d5ea5de9f96
# ╠═9b60de2c-19e6-40ca-b4b0-bc60b19f3adb
# ╟─c9f608c8-0088-4667-909e-2ed0cc0d1866
# ╟─78a5a17f-8064-437f-bdf5-6ac47fe0741c
