using Pkg; Pkg.activate()
using Infiltrator, Revise
using Pkg; Pkg.activate("./scripts")
using RootUntangling
using JuMP, HiGHS, Gurobi
using CairoMakie
using Dates, Statistics
import .Makie: HSV
ENV["JULIA_DEBUG"] = RootUntangling

# choose boy
playground_dir = "playground"
difficulty = "easy"
validation_dir = "validation/$(difficulty)"

directory = validation_dir
roi_nr = 1

# read data

begin
    today = now() |> monthday .|> string .|> (x -> length(x) == 1 ? "0" * x : x) |> x -> x[1] * "-" * x[2]

    cluster_threshold = 3
    loop_threshold = 5
    reverse_y = true
    angle_correction = -

    filename_segments = "./data/$(directory)/ROI_$(roi_nr)/segment_info_with_coords.csv"
    filename_vertices = "./data/$(directory)/ROI_$(roi_nr)/bp1_segments_grouped.csv"
    pg = get_pregraph(filename_segments, filename_vertices; cluster_threshold, loop_threshold, reverse_y, angle_correction)
    rg_full = get_rootgraph(pg)

    graphplot(rg_full, augmented_alpha = 0.01)
end

begin
    min_vertices = 20
    y_min = -1000
    y_max = 4800

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

subidx = 2
rg = rgs[subidx]
graphplot(rg, augmented_alpha = 0.01)

# NN predictions
lines(
    rgs[subidx], E₀(rgs[subidx]),
    color = [fill(HSV(0, 1, pred_primary(re)), 3) for re in E₀(rgs[subidx])] |> x -> reduce(vcat, x)
)
# savefig(homedir() * "\\Downloads\\wa.svg")

model, time = @timed solve_rsa(
    rg; optimizer = HiGHS.Optimizer, time_limit = 2*60,
    num_roots = 1, max_overlapping = 5,
    # ρₘ_max = 0.5, ρₘ_min = 0, ρᵧ_max = 0.5, ρᵧ_min = 0,
    # ρₛ = 0.5
)

f_g = graphplot(rg, model, augmented_alpha = 0.3, size = (600, 800), vertex_kwargs = Dict(:markersize => 2))
rss_raw = get_rootsystems(rg, pg, model);
rss = greedy_switch(rss_raw, f_obj = weighted_tortuosity);

r1 = rootplot(rss_raw, size = (600, 600), title = "Weighted tortuosity: $(weighted_tortuosity(rss_raw))")
r2 = rootplot(rss, size = (600, 600), title = "Weighted tortuosity: $(weighted_tortuosity(rss))")
# save(homedir() * "/Downloads/test1.svg", r1)
# save(homedir() * "/Downloads/test2.svg", r2)

# ap = annotation_plot(rg_full, size = (2000, 4000));
# save(homedir() * "/Downloads/thefuckingmimeryyysoto.svg", ap)
# ap = annotation_plot(rg, size = (2000, 4000));
# save(homedir() * "/Downloads/thefuckingmimeryyy.svg", ap)
# f_g = graphplot(rg, model, augmented_alpha = 0.3, size = (1200, 1500), vertex_kwargs = Dict(:markersize => 2))
# save(homedir() * "/Downloads/woagh.svg", f_g)
# r1 = rootplot(rss, size = (1200, 1500), title = "Weighted tortuosity: $(weighted_tortuosity(rss))")
# save(homedir() * "/Downloads/aghwuuu.svg", r1)




# validation

# categories = ["baby", "easy", "intermediate", "tough"]
# n_rois = [3, 6, 6, 4]
# for (i, category) in enumerate(categories), roi in 1:n_rois[i]
#     file = "data/validation/$(category)/ROI_$(roi)/GT.txt"
#     fix_groundtruth(file)
# end

filename_GT = "data/validation/$(difficulty)/ROI_$(roi_nr)/GT.txt"
vrss, segmentation_errors = get_validation_data(filename_GT, filename_segments, filename_vertices;
    cluster_threshold, loop_threshold, reverse_y, angle_correction
);
rootplot(vrss)

metric = standard_segments
rs_mapping, root_mappings = get_optimal_mappings(rss, vrss, metric, optimizer = HiGHS.Optimizer)

🔴(x) = round(Int64, x)
open("test.txt", "w") do file
    write(file, "class, roi, nr, type, PR_length, LR_count, LR_lengths")
end

figures = Makie.Figure[]
for (i, rs) in enumerate(rss)
    vrs = vrss[rs_mapping[i]]

    open("test.txt", "a") do file
        write(
            file,
            "\n$(difficulty), $(roi_nr), $i, GT, $(🔴(curve_length(primary(vrs)))), " *
                "$(🔴(length(laterals(vrs)))), " *
                "\"" * "$(🔴.(curve_length.(laterals(vrs))))"[2:end-1] * "\""
        )
    end

    open("test.txt", "a") do file
        write(
            file,
            "\n$(difficulty), $(roi_nr), $i, predicted, $(🔴(curve_length(primary(rs)))), " *
                "$(🔴(length(laterals(rs)))), " *
                "\"" * "$(🔴.(curve_length.(laterals(rs))))"[2:end-1] * "\""
        )
    end

    rs_copy, vrs_copy = deepcopy(rs), deepcopy(vrs);
    add_dummy_roots!(rs_copy, vrs_copy)
    final_loss = get_loss(rs_copy, vrs_copy, root_mappings[i], metric, verbose = true)
    laterals(vrs_copy) .= laterals(vrs_copy)[root_mappings[i]];

    f = rootplot(vrs_copy, rs_copy)
    Axis(f[1, 1]).title = "GT"
    Axis(f[1, 3]).title = "Predicted"
    Label(f[0, :], "loss = $(final_loss) $(string(metric))")
    push!(figures, f)
    save(homedir() * "/Downloads/validation_$(difficulty)_$(roi_nr)_$i.svg", f)
end

figures[1]
rss[1]
vrss[3]

rootplot(rs_copy)
rootplot(vrs_copy)

import RootUntangling: full_setdiff

[
    full_setdiff(standard_segments(laterals(rs_copy)[i]), standard_segments(laterals(vrs_copy)[i]))
    for i in eachindex(laterals(rs_copy))
]
