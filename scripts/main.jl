using Pkg; Pkg.activate()
using Infiltrator, Revise
using Pkg; Pkg.activate("./scripts")
using RootUntangling, RootUntangling.Plots
using JuMP
using HiGHS, Gurobi
using Dates
ENV["JULIA_DEBUG"] = RootUntangling

# choose boy
roi_nr = 12

# read data

begin
    today = now() |> monthday .|> string .|> (x -> length(x) == 1 ? "0" * x : x) |> x -> x[1] * "-" * x[2]

    nₕ_min = 1
    pₛ = 0.2
    dist_threshold = 3
    reverse_y = true

    filename_segments = "./data/ROI_$(roi_nr)/segment_info_with_coords.csv"
    filename_vertices = "./data/ROI_$(roi_nr)/bp1_segments_grouped.csv"
    rg = get_rootgraph(filename_segments, filename_vertices; dist_threshold, reverse_y, pₛ, nₕ_min)

    hypothesis_plot(rg)
end

model, time = @timed solve_rsa(
    rg; optimizer = Gurobi.Optimizer, add_momentum = true, time_limit = 10 * 60, hotstart_time = 2 * 60,
    num_roots = 2, ρₐ = 0.01, ρₕ = 0.97, ρₘ_max = 0.75, ρₙₙ_max = 0.9, ρᵧ_max = 0.5
)

roots = get_roots(rg, model);
plot(roots, size = (800, 800), title = "Time: $(round(time / 60, digits = 1)) min")
examine(roots)

savefig(homedir() * "/Downloads/test.svg")
savefig("results/roi$(roi_nr)_roots_$(today).svg")

# multi

rgs = get_subgraphs(rg; pₛ, nₕ_min) |>
    rgs -> filter(x -> length(V₀(x)) > 20, rgs);

begin
    plot(legend = false)
    for (i, rg) in enumerate(rgs)
        plot!(rg, color = HSV(i / length(rgs) * 360, 1, 1), augmented_alpha = 0.05)
    end
    plot!()
end

subidx = 1

plot(rgs[subidx], size = (1000, 800))
hypothesis_plot(rgs[subidx])

model, time = @timed solve_rsa(
    rgs[subidx]; optimizer = Gurobi.Optimizer, time_limit = 13 * 60, hotstart_time = 2 * 60,
    num_roots = 1
)

# plot(rgs[subidx], get_segment_classification_dict(rgs[subidx], model), size = (800, 800))
# savefig("results/roi$(roi_nr)-$(subidx)_classification_$(today).svg")

roots = get_roots(rgs[subidx], model)
plot(roots, size = (800, 800), title = "Time: $(round(time / 60, digits = 1)) min", lw = 1)
savefig("results/roi$(roi_nr)-$(subidx)_roots_$(today).svg")

# NN predictions

plot(
    rgs[subidx], size = (800, 800),
    edge_kwargs = Dict(:color => [HSV(0, 1, pred_primary(seg)) for seg in segments(rgs[subidx])] |> x -> reshape(x, 1, :)),
    vertex_kwargs = Dict(:color => [HSV(120, 1, pred_split(rv)) for rv in V(rgs[subidx])])
)
savefig(homedir() * "\\Downloads\\wa.svg")

# testing grounds

## tort tests
i = 1
begin
    r = roots[i]
    i += 1
    plot(r, title = "$(RootUntangling.tortuosity(r))")
end

## does switching work
isdefined(Main, :rgs) && (rg = rgs[subidx]);
begin
    roots = get_roots(rg, model)
    overlap_dict = RootUntangling.find_overlaps(rg, model, roots)
    overlap_segments = collect(keys(overlap_dict))
    seg = overlap_segments[1]
    r1 = overlap_dict[seg][1]
    r2 = overlap_dict[seg][2]

    p_before = plot(roots, size = (1000, 800), title = "Total tort: $(tortuosity(roots))")
    plot!(rg, seg, color = :red, linestyle = :solid, lw = 5, alpha = 0.3)
    RootUntangling.switch!(rg, seg, r1, r2)
    p_after = plot(roots, title = "Total tort: $(tortuosity(roots))")
    plot!(rg, seg, color = :red, linestyle = :solid, lw = 5, alpha = 0.3, label = false)
    plot(p_before, p_after)
end

## does greedy search work
isdefined(Main, :rgs) && (rg = rgs[subidx]);
begin
    roots = get_roots(rg, model)
    roots_improved = greedy_switch(rg, model, roots; max_tries = 100)

    p_before = plot(roots, size = (1000, 800), title = "Total tort: $(tortuosity(roots))")
    p_after = plot(roots_improved, title = "Total tort: $(tortuosity(roots_improved))")
    plot(p_before, p_after)
end
savefig(homedir() * "/Downloads/tortitup.svg")
