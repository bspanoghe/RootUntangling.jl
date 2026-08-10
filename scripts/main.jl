using Pkg; Pkg.activate()
using Infiltrator, Revise
using Pkg; Pkg.activate("./scripts")
using RootUntangling, RootUntangling.Plots
using JuMP
using HiGHS, Gurobi

using Dates
today = now() |> monthday .|> string .|> (x -> length(x) == 1 ? "0"*x : x) |> x -> x[1] * "-" * x[2]

# choose boy
roi_nr = 5

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

histogram(length.(vertices.(Vₕ₀(sg))))

model, time =  @timed RootUntangling.solve_rsa(
    sg; optimizer = Gurobi.Optimizer, add_momentum = true, time_limit = 13*60, hotstart_time = 2*60, 
    num_roots = 1, ρₐ = 0.01, ρₕ = 0.9
)

roots = get_roots(sg, model);
plot(roots, size = (800, 800), title = "Time: $(round(time/60, digits = 1)) min")
savefig(homedir() * "/Downloads/test.svg")
savefig("results/roi$(roi_nr)_roots_$(today).svg")

# multi

sgs = get_subgraphs(sg; pₛ, nₕ_min);
sgs = filter(x -> length(x) > 50, sgs);

begin
    plot(legend = false)
    for (i, sg) in enumerate(sgs)
        plot!(sg, color = HSV(i/length(sgs) * 360, 1, 1), augmented_alpha = 0.05)
    end
    plot!()
end

subidx = 2
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
plot(roots, size = (800, 800), title = "Time: $(round(time/60, digits = 1)) min", lw = 1)
savefig("results/roi$(roi_nr)-$(subidx)_roots_$(today).svg")

# NN predictions

plot(
    sgs[subidx], size = (800, 800),
    edge_kwargs = Dict(:color => [HSV(0, 1, pred_primary(he)) for he in Eₕ(sgs[subidx])] |> x -> reshape(x, 1, :)),
    vertex_kwargs = Dict(:color => [HSV(120, 1, pred_split(hv)) for hv in Vₕ(sgs[subidx])])
)
savefig(homedir() * "\\Downloads\\wa.svg")

# testing grounds

i = 1
begin
    r = roots[i]
    i += 1
    plot(r, title = "$(RootUntangling.tortuosity(r))")
end

he_classification_dict = RootUntangling.get_he_classification_dict(sg, model)
se_classification_dict = RootUntangling.get_se_classification_dict(sg, model)

edge_groups = [
    [
        se
        for se in E(sg, he)
        if abs(se_classification_dict[se]) > 0
    ] 
    for he in Eₕ₀(sg)
] |> v -> filter(x -> length(x) >= 2, v)

begin
    rr = roots[3]
    plot(roots, label = nothing)
    plot!(sg, reduce(vcat, edge_groups), color = :pink, lw = 5, alpha = 0.3, label = nothing)
    annotate!([(xs(rr)[i], ys(rr)[i], "$i", 6) for i in eachindex(xs(rr))])
end