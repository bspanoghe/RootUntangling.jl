using Pkg; Pkg.activate()
using Infiltrator, Revise
using Pkg; Pkg.activate("./scripts")
using RootUntangling
using JuMP, HiGHS, Gurobi
using CairoMakie
using Dates, Statistics
ENV["JULIA_DEBUG"] = RootUntangling

# choose boy

playground_dir = "playground"
difficulty = "tough"
validation_dir = "validation/$(difficulty)"

directory = playground_dir
roi_nr = 1

# read data

begin
    today = now() |> monthday .|> string .|> (x -> length(x) == 1 ? "0" * x : x) |> x -> x[1] * "-" * x[2]

    dist_threshold = 3
    reverse_y = true

    filename_segments = "./data/$(directory)/ROI_$(roi_nr)/segment_info_with_coords.csv"
    filename_vertices = "./data/$(directory)/ROI_$(roi_nr)/bp1_segments_grouped.csv"
    rg_full = get_rootgraph(filename_segments, filename_vertices; dist_threshold, reverse_y)

    graphplot(rg_full, augmented_alpha = 0.01)
end

begin
    min_vertices = 10
    y_threshold = -1000

    rgs = get_subgraphs(rg_full) |>
        rgs -> filter(rg -> length(rg) > min_vertices, rgs) |>
        rgs -> filter(rg -> minimum(y.(V₀(rg))) > y_threshold, rgs) |>
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

subidx = 1
rg = rgs[subidx]
graphplot(rg)

begin
    model, time = @timed solve_rsa(
        rg; optimizer = HiGHS.Optimizer, time_limit = 15*60,
        num_roots = 1, max_roots = 5
    )

    annotate_that_thang = false
    if annotate_that_thang
        f_g = graphplot(rg, model, augmented_alpha = 0.3, size = (1000, 2000))
        rd = get_result_dict(rg, model)
        annotation_coords = [(mean(xs(rg, re)), mean(ys(rg, re))) for re in E₀(rg)]
        annotation_texts = [
            "($(rd[re][:e₊]) / $(rd[re][:e₋]))"
            for re in E₀(rg)
        ]
        annotation!(f_g.content[1], annotation_coords; text = annotation_texts, color = :red, shrink = (0, 0), fontsize = 8)

        annotation_coords = [(x(rv), y(rv)) for rv in V₀(rg)]
        annotation_texts = string.(id.(V₀(rg)))
        annotation!(f_g.content[1], annotation_coords; text = annotation_texts, color = :green, shrink = (0, 0), fontsize = 8)
        save(homedir() * "/Downloads/oooo_the_mimeryyy.svg", f_g)

    else
       f_g = graphplot(rg, model, augmented_alpha = 0.3, size = (200, 500), vertex_kwargs = Dict(:markersize => 2))
       f_g
    end
end
rss = get_rootsystems(rg, model);
rss_new = greedy_switch(rg, rss, f_obj = weighted_tortuosity);

r1 = rootplot(rg, rss, size = (600, 600), title = "Weighted tortuosity: $(weighted_tortuosity(rg, rss))")
r2 = rootplot(rg, rss_new, size = (600, 600), title = "Weighted tortuosity: $(weighted_tortuosity(rg, rss_new))")

save(homedir() * "/Downloads/test1.svg", r1)
save(homedir() * "/Downloads/test2.svg", r2)

# NN predictions
import .Makie: HSV
lines(
    rgs[subidx], E₀(rgs[subidx]),
    color = [fill(HSV(0, 1, pred_primary(re)), 3) for re in E₀(rgs[subidx])] |> x -> reduce(vcat, x)
)
savefig(homedir() * "\\Downloads\\wa.svg")