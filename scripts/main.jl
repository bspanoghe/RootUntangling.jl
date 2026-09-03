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
difficulty = "baby"
validation_dir = "validation/$(difficulty)"

directory = validation_dir
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
    sg_full = get_supergraph(filename_segments, filename_vertices; dist_threshold, reverse_y, pₛ, nₕ_min)

    hypothesis_plot(sg_full)
end

begin
    min_vertices = 5
    y_threshold = 1000

    sgs = get_subgraphs(sg_full; pₛ, nₕ_min) |>
        sgs -> filter(sg -> length(sg) > min_vertices, sgs) |>
        sgs -> filter(sg -> minimum(y.(Vₕ₀(sg))) > y_threshold, sgs) |>
        sgs -> sort(sgs, by = sg -> mean(x.(Vₕ₀(sg))));

    f_multi = Figure()
    ax_multi = Axis(f_multi[1, 1]; aspect = DataAspect())
    
    ls = [
        lines!(ax_multi, sg, Eₕ₀(sg), color = Makie.HSV(i / length(sgs) * 360, 1, 1))
        for (i, sg) in enumerate(sgs)
    ]
    Legend(f_multi[1, 2], ls, string.(1:length(ls)))

    f_multi
end

annotating = true
if annotating
    root_systems = get_root_systems(sgs, ones(Int64, length(sgs)), optimizer = Gurobi.Optimizer,
        time_limit = 13*60, hotstart_time = 2*60, ρₒ_base = 4.0
    )

    f_ann = annotation_plot(sg_full, root_systems, size = (3200, 1800))
    save("results/$(directory)/ROI_$(roi_nr).svg", f_ann)
    write_annotation("results/$(directory)/ROI_$(roi_nr).txt", sg_full, sgs, root_systems)
end

f_ann = annotation_plot(sg_full, root_systems, size = (1600, 900))

subidx = 1
sg = sgs[subidx]
hypothesis_plot(sg)

begin
    model, time = @timed solve_rsa(
        sg; optimizer = Gurobi.Optimizer, add_momentum = true, time_limit = 60, hotstart_time = 60,
        num_roots = 1, ρₒ_base = 4.0
    )

    roots = get_roots(sg, model);
    roots_new = greedy_switch(sg, model, roots)
end

r1 = rootplot(roots, size = (600, 600), title = "Time: $(round(time / 60, digits = 1)) min")
r2 = rootplot(roots_new, size = (600, 600), title = "Time: $(round(time / 60, digits = 1)) min")

save(homedir() * "/Downloads/test1.svg", r1)
save(homedir() * "/Downloads/test2.svg", r2)

# NN predictions

plot(
    sgs[subidx], size = (800, 800),
    edge_kwargs = Dict(:color => [HSV(0, 1, pred_primary(he)) for he in Eₕ(sgs[subidx])] |> x -> reshape(x, 1, :)),
    vertex_kwargs = Dict(:color => [HSV(120, 1, pred_split(hv)) for hv in Vₕ(sgs[subidx])])
)
savefig(homedir() * "\\Downloads\\wa.svg")

# testing grounds

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