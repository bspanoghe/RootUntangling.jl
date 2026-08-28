using Pkg; Pkg.activate()
using Infiltrator, Revise
using Pkg; Pkg.activate("./scripts")
using RootUntangling
using JuMP, HiGHS, Gurobi
using CairoMakie
using Dates, Statistics
ENV["JULIA_DEBUG"] = RootUntangling

# choose boy
directory = "validation/baby"
roi_nr = 1

# read data

begin
    today = now() |> monthday .|> string .|> (x -> length(x) == 1 ? "0" * x : x) |> x -> x[1] * "-" * x[2]

    nₕ_min = 1
    pₛ = 0.2
    dist_threshold = 3
    reverse_y = true

    filename_segments = "./data/$(directory)/ROI_$(roi_nr)/segment_info_with_coords.csv"
    filename_vertices = "./data/$(directory)/ROI_$(roi_nr)/bp1_segments_grouped.csv"
    sg = get_supergraph(filename_segments, filename_vertices; dist_threshold, reverse_y, pₛ, nₕ_min)

    hypothesis_plot(sg)
end

begin
    sgs = get_subgraphs(sg; pₛ, nₕ_min) |>
        sgs -> filter(x -> length(x) > 15, sgs) |>
        sgs -> sort(sgs, by = sg -> mean(x.(Vₕ₀(sg))));

    f_multi = Figure()
    ax_multi = Axis(f_multi[1, 1]; aspect = DataAspect())
    for (i, sg) in enumerate(sgs)
        lines!(ax_multi, sg, Eₕ₀(sg), color = Makie.HSV(i / length(sgs) * 360, 1, 1))
    end
    f_multi
end

for subidx in eachindex(sgs)
    sg = sgs[subidx]
    f_ann = annotation_plot(sg, size = (2400, 2000))
    save("results/$(directory)/ROI_$(roi_nr)_$(subidx).svg", f_ann)
end

model, time = @timed solve_rsa(
    sg; optimizer = Gurobi.Optimizer, add_momentum = true, time_limit = 60, hotstart_time = 60,
    num_roots = 1, ρₐ = 0.01, ρₘ_max = 0.75, ρₙₙ_max = 0.9, ρᵧ_max = 0.5, ρₒ_base = 4.0
)

roots = get_roots(sg, model);
graphplot(sg, get_he_classification_dict(sg, model))
r = rootplot(roots, height = 800, width = 600, title = "Time: $(round(time / 60, digits = 1)) min")
examine(roots)

save(homedir() * "/Downloads/test.png", r)
save("results/roi$(roi_nr)_roots_$(today).png", r)

annotation_dict = RootUntangling.get_annotation_dict(sg, roots)
RootUntangling.annotation_plot(sg, annotation_dict)

plot(roots, size = (800, 800), title = "Time: $(round(time / 60, digits = 1)) min", lw = 1)
savefig("results/roi$(roi_nr)-$(subidx)_roots_$(today).svg")

examine(roots)

# NN predictions

plot(
    sgs[subidx], size = (800, 800),
    edge_kwargs = Dict(:color => [HSV(0, 1, pred_primary(he)) for he in Eₕ(sgs[subidx])] |> x -> reshape(x, 1, :)),
    vertex_kwargs = Dict(:color => [HSV(120, 1, pred_split(hv)) for hv in Vₕ(sgs[subidx])])
)
savefig(homedir() * "\\Downloads\\wa.svg")

# testing grounds


## does greedy search work
isdefined(Main, :sgs) && (sg = sgs[subidx]);
begin
    roots = get_roots(sg, model)
    roots_improved = greedy_switch(sg, model, roots; max_tries = 100)

    p_before = plot(roots, size = (1000, 800), title = "Total roughness: $(roughness(roots))")
    p_after = plot(roots_improved, title = "Total roughness: $(roughness(roots_improved))")
    plot(p_before, p_after)
end
savefig(homedir() * "/Downloads/tortitup.svg")





points = Observable(Point2f[])

scene = Scene(camera = campixel!)
linesegments!(scene, points, color = :black)
scatter!(scene, points, color = :gray)

on(events(scene).mousebutton) do event
    if event.button == Mouse.left
        if event.action == Mouse.press || event.action == Mouse.release
            mp = events(scene).mouseposition[]
            push!(points[], mp)
            notify(points)
        end
    end
end

scene

working_on = 1

begin
    f = RootUntangling.annotation_plot(sg, annotation_dict)
    ax = Axis(f[1, 1])
    hidedecorations!(ax, label = false, ticks = false, ticklabels = false)
    hidespines!(ax, :t, :r)

    i = Observable((0.0,0.0))
    str = lift(i -> "$(i)", i)
    text!(ax, 1, -0.5, text = str,  align = (:center, :center))
    on(events(f).mousebutton, priority = 2) do event
        if event.button == Mouse.left && event.action == Mouse.press
            global elements = Makie.pick_sorted(f.scene, events(f).mouseposition[], 30)
            filter!(x -> x[1] isa Lines, elements)
            line_idx = findfirst(x -> x[1] isa Lines, elements)
            if !isnothing(line_idx)
                line_element = elements[line_idx]
                he_picked = line_element[1].arg2.value[]

                current_annotation_idx = findfirst(x -> x == working_on, annotation_dict[he_picked])
                if isnothing(current_annotation_idx)
                    push!(annotation_dict[he_picked], working_on)
                else
                    deleteat!(annotation_dict[he_picked], current_annotation_idx)
                end

                f = RootUntangling.annotation_plot(sg, annotation_dict)
            end
        end
    end
    f
end

f = RootUntangling.annotation_plot(sg, annotation_dict)
DataInspector(f)
f