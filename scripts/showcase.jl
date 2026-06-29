### A Pluto.jl notebook ###
# v0.20.21

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
	Pkg.activate("..")
	using RootUntangling
	using RootUntangling.Plots
end

# ╔═╡ 1423c220-fb6d-4c4c-bd07-140e29cc79cc
md"## Imports"

# ╔═╡ 71e10643-4c1d-4069-bd9e-eb8fd58dccc7
md"## Data reading"

# ╔═╡ 72b6f5ea-f063-4705-bc9c-b06d55fd0d43
num_hypotheses = 2

# ╔═╡ 02795186-93cf-40cb-af20-6869fa45ebdb
dist_threshold = 3

# ╔═╡ 3fc9f8b8-e6e8-4409-a920-98d30ac9c2d1
md"## Data visualisation"

# ╔═╡ 690ea173-af44-475e-bc5a-1f1044e48309
@bind roi_nr Select(1:7, default = 4)

# ╔═╡ c26c9948-8ff1-44c3-8b68-8e0907fdb3ec
begin
    filename_segments = "../clouddata/branchpoints/ROI_$(roi_nr)/segment_info_with_coords.csv";
    filename_vertices = "../clouddata/branchpoints/ROI_$(roi_nr)/bp1_segments_grouped.csv";
    sg = get_supergraph(filename_segments, filename_vertices; dist_threshold, num_hypotheses);
end

# ╔═╡ 0d47753a-4608-457f-a071-5e4e5d24d2ed
begin
	img_url = "../clouddata/branchpoints/ROI_$(roi_nr)/Skeleton.png"
	md"""$(LocalResource(img_url, :height => 500))"""
end

# ╔═╡ 498a192d-732e-4205-9043-402d5b4f3988
begin
	mrg = 50 # margin
	x_ll = 0-mrg
	x_ul = maximum(x.(V₀(sg)))+mrg
	y_ll = minimum(y.(V₀(sg)))-mrg
	y_ul = maximum(y.(V₀(sg)))+mrg
end;

# ╔═╡ c3115132-0524-40e7-87fd-e27a35375df7
@bind xmin Slider(x_ll:x_ul, show_value = true, default = x_ll)

# ╔═╡ 2adc1870-aad6-49e4-a38c-c13009111c48
@bind xmax Slider(x_ll:x_ul, show_value = true, default = x_ul)

# ╔═╡ b2655ff8-aa2e-4db9-ba2d-6e6a9463cdb3
@bind ymin Slider(y_ll:y_ul, show_value = true, default = y_ll)

# ╔═╡ e08b968e-f98e-4a48-9560-0cd01279ac38
@bind ymax Slider(y_ll:y_ul, show_value = true, default = y_ul)

# ╔═╡ 3f40ee53-68ad-4bfa-a784-bbad4a366f02
plot(sg,
	 legend = false, yflip = true, augmented_alpha = 0.0, aspect_ratio = :equal,
	 xlims = (xmin, xmax), ylims = (ymin, ymax)
)

# ╔═╡ 7894d17f-b98c-4e19-8d3b-4fb690aaeb42
md"## Solving"

# ╔═╡ 11335ca1-7f5b-4918-a45d-c28435315a73
@bind solver Select([solve_rsa, solve_rsa_simple, solve_rsa_fulloption])

# ╔═╡ 5355a78b-3bee-454e-a18b-0e6ac2d202b7
ep2f, el2f, vp2f, vl2f = solver(sg, time_limit = 120, α = 1, β = 2, γ = 0.8, δ = 2);

# ╔═╡ 3895237d-8864-4681-95d7-fd23563f23a5
#=╠═╡
begin
	classification_dict = [
		e => round(Bool, ep2f[e]) + round(Bool, el2f[e])*im 
		for e in E(sg)
	] |> Dict
	he_class_dict = [
		he => sum([classification_dict[e] for e in E(sg, he)])
		for he in Eₕ(sg)
	] |> Dict

	plot(
		sg, legend = false, yflip = true, augmented_alpha = 0.3, aspect_ratio = :equal, xlims = (xmin, xmax), ylims = (ymin, ymax),
		edge_kwargs = Dict(
			:color => [(real(he_class_dict[he]) > 0) * RGB(1.0, 0, 0) + (imag(he_class_dict[he]) > 0) * RGB(0, 0, 1.0) for he in Eₕ(sg)] |> x -> reshape(x, 1, :),
			:linewidth => [is_augmented(he) ? abs(he_class_dict[he]) : abs(he_class_dict[he]) + 0.5 for he in Eₕ(sg)] |> x -> reshape(x, 1, :)
		),
		vertex_kwargs = Dict(
			:color => :grey, :markersize => 2
		)
	)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─1423c220-fb6d-4c4c-bd07-140e29cc79cc
# ╠═bace35f3-a1e4-49ff-9f63-fb4ca65d2176
# ╠═7d5639a6-a788-4877-8ac8-394ffb45a220
# ╠═1ebe1d3c-dd67-44d5-b292-c6ce4b9f9aa1
# ╟─71e10643-4c1d-4069-bd9e-eb8fd58dccc7
# ╠═72b6f5ea-f063-4705-bc9c-b06d55fd0d43
# ╠═02795186-93cf-40cb-af20-6869fa45ebdb
# ╠═c26c9948-8ff1-44c3-8b68-8e0907fdb3ec
# ╟─3fc9f8b8-e6e8-4409-a920-98d30ac9c2d1
# ╠═690ea173-af44-475e-bc5a-1f1044e48309
# ╟─0d47753a-4608-457f-a071-5e4e5d24d2ed
# ╟─3f40ee53-68ad-4bfa-a784-bbad4a366f02
# ╟─498a192d-732e-4205-9043-402d5b4f3988
# ╠═c3115132-0524-40e7-87fd-e27a35375df7
# ╠═2adc1870-aad6-49e4-a38c-c13009111c48
# ╠═b2655ff8-aa2e-4db9-ba2d-6e6a9463cdb3
# ╠═e08b968e-f98e-4a48-9560-0cd01279ac38
# ╟─7894d17f-b98c-4e19-8d3b-4fb690aaeb42
# ╠═11335ca1-7f5b-4918-a45d-c28435315a73
# ╠═5355a78b-3bee-454e-a18b-0e6ac2d202b7
# ╟─3895237d-8864-4681-95d7-fd23563f23a5
