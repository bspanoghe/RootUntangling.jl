using Pkg; Pkg.activate()
using Infiltrator, Revise
using Pkg; Pkg.activate("./scripts")
using RootUntangling, RootUntangling.Plots
using JuMP
using HiGHS, Gurobi

using Dates
today = now() |> monthday .|> string .|> (x -> length(x) == 1 ? "0"*x : x) |> x -> x[1] * "-" * x[2]

# choose boy
roi_nr = 15

# read data

begin
    nₕ_min = 1
    pₛ = 0.2
    dist_threshold = 3
    flip_y = true

    filename_segments = "./data/ROI_$(roi_nr)/segment_info_with_coords.csv";
    filename_vertices = "./data/ROI_$(roi_nr)/bp1_segments_grouped.csv";
    sg = get_supergraph(filename_segments, filename_vertices; dist_threshold, flip_y, pₛ, nₕ_min);

    hypothesis_plot(sg)
end

length(V₀(sg))
histogram(length.(vertices.(Vₕ₀(sg))))

model, time =  @timed RootUntangling.solve_rsa(
    sg; optimizer = Gurobi.Optimizer, add_momentum = true, time_limit = 13*60, hotstart_time = 2*60, 
    num_roots = 1, ρₐ = 0.01, ρₕ = 0.9
)

roots = get_roots(sg, model);
plot(roots, size = (800, 800), title = "Time: $(round(time/60, digits = 1)) min")
savefig("results/roi$(roi_nr)_roots_$(today).svg")

# multi

sgs = get_subgraphs(sg; pₛ, nₕ_min)
sgs = filter(x -> length(x) > 50, sgs)

begin
    plot(legend = false)
    for (i, sg) in enumerate(sgs)
        plot!(sg, color = HSV(i/length(sgs) * 360, 1, 1), augmented_alpha = 0.05)
    end
    plot!()
end

subidx = 1
plot(sgs[subidx], size = (1000, 800))
savefig(homedir() * "\\Downloads\\wa.svg")
hypothesis_plot(sgs[subidx])

model, time = @timed RootUntangling.solve_rsa(
    sgs[subidx]; optimizer = Gurobi.Optimizer, add_momentum = true, time_limit = 13*60, hotstart_time = 2*60,
    num_roots = 1, ρₐ = 0.01, ρₕ = 0.9, ρₘ_max = 0.9
)

plot(sgs[subidx], get_he_classification_dict(sgs[subidx], model), size = (800, 800))
savefig("results/roi$(roi_nr)-$(subidx)_classification_$(today).svg")

roots = get_roots(sgs[subidx], model)
plot(roots, size = (800, 800), title = "Time: $(round(time/60, digits = 1)) min", linewidth = 0.5)
savefig("results/roi$(roi_nr)-$(subidx)_roots_$(today).svg")

# NN predictions

plot(
    sgs[subidx], size = (800, 800),
    edge_kwargs = Dict(:color => [HSV(0, 1, pred_primary(he)) for he in Eₕ(sgs[subidx])] |> x -> reshape(x, 1, :)),
    vertex_kwargs = Dict(:color => [HSV(120, 1, pred_split(hv)) for hv in Vₕ(sgs[subidx])])
)
savefig(homedir() * "\\Downloads\\wa.svg")
