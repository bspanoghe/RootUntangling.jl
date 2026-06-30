using Pkg; Pkg.activate()
using Infiltrator, Revise
using Pkg; Pkg.activate("./scripts")
using RootUntangling, RootUntangling.Plots
using JuMP
using HiGHS, Gurobi

using Dates
time = now() |> monthday .|> string .|> (x -> length(x) == 1 ? "0"*x : x) |> x -> x[1] * "-" * x[2]

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

    nₕs = [
        maximum([length(vertices(hv)) for hv in RootUntangling.gethypervertex.(vertices(he), [Vₕ(sg)])])
        for he in Eₕ(sg)
    ]
    nₕ₀s = [
        maximum([length(vertices(hv)) for hv in RootUntangling.gethypervertex.(vertices(he), [Vₕ(sg)])])
        for he in Eₕ₀(sg)
    ]

    ΔH = 80

    plot(sg, size = (600, 800), label = false, edge_kwargs = 
        Dict(
            :linewidth => [width(he)/2 for he in Eₕ(sg)] |> x -> reshape(x, 1, :),
            :linecolor => [HSV(ΔH*nₕ, 1, 0.75) for nₕ in nₕs] |> x -> reshape(x, 1, :),
        ),
    )
    plot!(
        fill(missing, 1, length(unique(nₕs))), fill(missing, 1, length(unique(nₕs))), 
        lw = 2, legend = true, legendfontsize = 12,
        label = reshape(sort(unique(nₕs)), 1, :), 
        linecolor = reshape([HSV(ΔH*nₕ, 1, 0.75) for nₕ in sort(unique(nₕs))], 1, :)
    )
    plot!(sg, [he for (nₕ₀, he) in zip(nₕ₀s, Eₕ₀(sg)) if nₕ₀ == 1], aspect_ratio = :equal,
        label = false, color = HSV(ΔH, 1, 0.75), lw = 5
    )
end

length(V₀(sg))
histogram(nₕ₀s)

@time model = RootUntangling.solve_rsa(
    sg; optimizer = Gurobi.Optimizer, add_momentum = true, time_limit = 13*60, hotstart_time = 2*60, 
    num_roots = 1, ρₐ = 0.01, ρₕ = 0.9
)

plot(sg, get_he_classification_dict(sg, model), augmented_alpha = 0.1, size = (800, 800))
savefig("results/roi$(roi_nr)_classification_$(time).svg")

roots = get_roots(sg, model);
plot(roots, size = (800, 800))
savefig("results/roi$(roi_nr)_roots_$(time).svg")

# multi

sgs = get_subgraphs(sg)
sgs = filter(x -> length(x) > 30, sgs)

begin
    plot(legend = false)
    for (i, sg) in enumerate(sgs)
        plot!(sg, color = HSV(i/length(sgs) * 360, 1, 1), augmented_alpha = 0.05)
    end
    plot!()
end

subidx = 1
plot(sgs[subidx], size = (1000, 800))

model = RootUntangling.solve_rsa(
    sgs[subidx]; optimizer = Gurobi.Optimizer, add_momentum = true, time_limit = 10*60, hotstart_time = 1*60,
    num_roots = 2, ρₐ = 0.01, ρₕ = 0.9, ρₘ_max = 0.75, ρₙₙ_max = 0.99
)

plot(sgs[subidx], get_he_classification_dict(sgs[subidx], model), size = (800, 800))
savefig("results/roi$(roi_nr)-$(subidx)_classification_$(time).svg")

roots = get_roots(sgs[subidx], model)
plot(roots, size = (800, 800))
savefig("results/roi$(roi_nr)-$(subidx)_roots_$(time).svg")

# NN predictions

[pred_primary(he) for he in Eₕ₀(sg)]

plot(
    sgs[subidx], size = (800, 800),
    edge_kwargs = Dict(:color => [HSV(0, 1, pred_primary(he)) for he in Eₕ(sgs[subidx])] |> x -> reshape(x, 1, :)),
    vertex_kwargs = Dict(:color => [HSV(120, 1, pred_split(hv)) for hv in Vₕ(sgs[subidx])])
)

# angle tests

begin
    v = V₀(sg)[end-6]
    plot(aspect_ratio = :equal)
    plot!(v)
    for he in Eₕ₀(v)
        e = E(sg, he)[1]
        println(vertices(e))
        α = RootUntangling.angle(sg, e, reverse_order = false)
        # α = RootUntangling.cosine_similarity(sg, e, -pi/2, reverse_order = true)
        plot!(sg, e, label = "$(round(α, digits = 2))")
    end
    plot!()
end

# width tests

Eₕ₀(sg) .|> width |> x -> quantile(x, 0.9)
