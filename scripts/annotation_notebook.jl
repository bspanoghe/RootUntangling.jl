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
	using PlutoUI; TableOfContents()
end

# ╔═╡ 1ebe1d3c-dd67-44d5-b292-c6ce4b9f9aa1
begin
	Pkg.activate(".")
	using RootUntangling
	using WGLMakie
	using JuMP, Gurobi
end

# ╔═╡ 1423c220-fb6d-4c4c-bd07-140e29cc79cc
md"## Imports"

# ╔═╡ 71e10643-4c1d-4069-bd9e-eb8fd58dccc7
md"## Data reading"

# ╔═╡ 72b6f5ea-f063-4705-bc9c-b06d55fd0d43
begin
    nₕ_min = 1
    pₛ = 0.2
    dist_threshold = 3
    reverse_y = true
end;

# ╔═╡ 3fc9f8b8-e6e8-4409-a920-98d30ac9c2d1
md"## Data visualisation"

# ╔═╡ 690ea173-af44-475e-bc5a-1f1044e48309
@bind roi_nr Select(1:7, default = 4)

# ╔═╡ c26c9948-8ff1-44c3-8b68-8e0907fdb3ec
begin
    filename_segments = "../data/ROI_$(roi_nr)/segment_info_with_coords.csv"
    filename_vertices = "../data/ROI_$(roi_nr)/bp1_segments_grouped.csv"
    rg = get_rootgraph(filename_segments, filename_vertices; dist_threshold, reverse_y, pₛ, nₕ_min)
end

# ╔═╡ 3f40ee53-68ad-4bfa-a784-bbad4a366f02
hypothesis_plot(rg)

# ╔═╡ 7894d17f-b98c-4e19-8d3b-4fb690aaeb42
md"## Solving"

# ╔═╡ 990ad7fb-4256-41ce-a6fe-520a103c5ea7
begin
	time_limit = 60
	hotstart_time = 60
	num_roots = 1
end;

# ╔═╡ 5355a78b-3bee-454e-a18b-0e6ac2d202b7
# ╠═╡ show_logs = false
model, time = @timed solve_rsa(
    rg; optimizer = Gurobi.Optimizer, time_limit, hotstart_time, num_roots
);

# ╔═╡ 797a5241-5c5c-47f9-a79c-c1ceec3add21
roots = get_roots(rg, model);

# ╔═╡ 3ab99e7a-2e94-486c-8f25-6ca8aaad21b9
r = rootplot(roots, title = "Time: $(round(time / 60, digits = 1)) min")

# ╔═╡ 95fce8b3-1db5-4a61-b4e3-0e8e36181416
annotation_dict = RootUntangling.get_annotation_dict(rg, roots);

# ╔═╡ 08de3104-630e-4792-836f-0ab11d324470
@bind working_on Select(
	sort(unique(reduce(vcat, values(annotation_dict))))
)

# ╔═╡ d9658eec-968d-44d1-bbb7-21717e4a4bfd
@bind clicked CounterButton("Update")

# ╔═╡ 2a4561ca-f05e-4fba-b14a-2d5ea5de9f96
begin
    clicked
    f = RootUntangling.annotation_plot(rg, annotation_dict)
    ax = Axis(f[1, 1])

    if !@isdefined annotation_dict_updated
        annotation_dict_updated = annotation_dict
    end

    on(events(f).mousebutton, priority = 2) do event
        if event.button == Mouse.left && event.action == Mouse.press
            elements = Makie.pick_sorted(f.scene, events(f).mouseposition[], 30)
            filter!(x -> x[1] isa Lines, elements)
            line_idx = findfirst(x -> x[1] isa Lines, elements)
            if !isnothing(line_idx)
                line_element = elements[line_idx]
                re_picked = line_element[1].arg2.value[]

                current_annotation_idx = findfirst(
                    x -> x == working_on, annotation_dict[re_picked]
                )
                if isnothing(current_annotation_idx)
                    annotation_dict_updated = push!(annotation_dict[re_picked], working_on)
                else
                    annotation_dict_updated = deleteat!(annotation_dict[re_picked], current_annotation_idx)
                end
            end
        end
    end
    f
end

# ╔═╡ Cell order:
# ╟─1423c220-fb6d-4c4c-bd07-140e29cc79cc
# ╠═bace35f3-a1e4-49ff-9f63-fb4ca65d2176
# ╠═7d5639a6-a788-4877-8ac8-394ffb45a220
# ╠═1ebe1d3c-dd67-44d5-b292-c6ce4b9f9aa1
# ╟─71e10643-4c1d-4069-bd9e-eb8fd58dccc7
# ╠═72b6f5ea-f063-4705-bc9c-b06d55fd0d43
# ╠═c26c9948-8ff1-44c3-8b68-8e0907fdb3ec
# ╟─3fc9f8b8-e6e8-4409-a920-98d30ac9c2d1
# ╠═690ea173-af44-475e-bc5a-1f1044e48309
# ╠═3f40ee53-68ad-4bfa-a784-bbad4a366f02
# ╟─7894d17f-b98c-4e19-8d3b-4fb690aaeb42
# ╠═990ad7fb-4256-41ce-a6fe-520a103c5ea7
# ╠═5355a78b-3bee-454e-a18b-0e6ac2d202b7
# ╠═797a5241-5c5c-47f9-a79c-c1ceec3add21
# ╠═3ab99e7a-2e94-486c-8f25-6ca8aaad21b9
# ╠═95fce8b3-1db5-4a61-b4e3-0e8e36181416
# ╟─08de3104-630e-4792-836f-0ab11d324470
# ╟─d9658eec-968d-44d1-bbb7-21717e4a4bfd
# ╟─2a4561ca-f05e-4fba-b14a-2d5ea5de9f96
