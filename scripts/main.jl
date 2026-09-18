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

directory = playground_dir
roi_nr = 6

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

annotating = false
if annotating
    root_systems = get_root_systems(rgs, ones(Int64, length(rgs)), optimizer = Gurobi.Optimizer,
        time_limit = 13*60, hotstart_time = 2*60, ρₒ_base = 4.0
    )

    f_ann = annotation_plot(rg_full, root_systems, size = (3200, 1800))
    save("results/$(directory)/ROI_$(roi_nr).svg", f_ann)
    write_annotation("results/$(directory)/ROI_$(roi_nr).txt", rg_full, rgs, root_systems)
end

# f_ann = annotation_plot(rg_full, root_systems, size = (1600, 900))

subidx = 1
rg = rgs[subidx]
graphplot(rg)

begin
    model, time = @timed solve_rsa(
        rg; optimizer = HiGHS.Optimizer, time_limit = 60,
        num_roots = 1, ρₘ_max = 0.5
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
    else
        f_g = graphplot(rg, model, augmented_alpha = 0.3, size = (250, 500), vertex_kwargs = Dict(:markersize => 2))
    end

    f_g
end
save(homedir() * "/Downloads/wwawa.svg", f_g)

roots = get_rootsystems(rg, model);
roots_new = greedy_switch(rg, model, roots)


r1 = rootplot(roots, size = (600, 600), title = "Time: $(round(time / 60, digits = 1)) min")
r2 = rootplot(roots_new, size = (600, 600), title = "Time: $(round(time / 60, digits = 1)) min")

save(homedir() * "/Downloads/test1.svg", r1)
save(homedir() * "/Downloads/test2.svg", r2)

# NN predictions
import .Makie: HSV
lines(
    rgs[subidx], E₀(rgs[subidx]),
    color = [fill(HSV(0, 1, pred_primary(re)), 3) for re in E₀(rgs[subidx])] |> x -> reduce(vcat, x)
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
    f = RootUntangling.annotation_plot(rg, annotation_dict)
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
                re_picked = line_element[1].arg2.value[]

                current_annotation_idx = findfirst(x -> x == working_on, annotation_dict[re_picked])
                if isnothing(current_annotation_idx)
                    push!(annotation_dict[re_picked], working_on)
                else
                    deleteat!(annotation_dict[re_picked], current_annotation_idx)
                end

                f = RootUntangling.annotation_plot(rg, annotation_dict)
            end
        end
    end
    f
end

f = RootUntangling.annotation_plot(rg, annotation_dict)
DataInspector(f)
f




###########################################################################################################
T, U = typeof(rg).parameters

result_dict = get_result_dict(rg, model)

# get all primary root fragments
primaries = Root{T, U}[]

c_counts = Dict(E₂(rg) .=> round.(Int64, value.(model[:cp])))
e₊_counts = Dict(E(rg) .=> round.(Int64, value.(model[:ep₊])))
e₋_counts = Dict(E(rg) .=> round.(Int64, value.(model[:ep]) - value.(model[:ep₊])))
# filter!(x -> x.second > 0, c_counts) #!

while !isempty(c_counts)
    root = grow_root!(c_counts, :primary)
end

# assign laterals that split to their primary root

# assign laterals that appeared to their most probable primary root




# get primary root(s)




import RootUntangling: RootEdge

function grow_root!(c_counts::Dict, e₊_counts::Dict, e₋_counts::Dict, rg::RootGraph, root_type::Symbol)
    root_type ∈ [:primary, :lateral] || error("Root type should be primary or lateral")
    
    # start with a random active connection
    c0 = findfirst(x -> x > 0, c_counts)
    c_counts[c0] -= 1

    # determine what directions are possible for the two edges
    es = c0
    follows_polarity = Bool[]
    opposite_direction = allequal(src.(es)) || allequal(dst.(es)) # both edges point towards or away from shared vertex
    if opposite_direction # edges require different polarity
        if e₊_counts[es[1]] > 0 && e₋_counts[es[2]] > 0
            append!(follows_polarity, [true, false])
        else
            append!(follows_polarity, [false, true])
        end
    else # edges require same polarity
        if e₊_counts[es[1]] > 0 && e₊_counts[es[2]] > 0
            append!(follows_polarity, [true, true])
        else
            append!(follows_polarity, [false, false])
        end
    end
    
    has_restarted = false
    fullgrown = false
    # until you can no longer append root segments in either direction:
    while !fullgrown
        # choose an edge at an end of the root
        e_current = has_restarted ? es[1] : es[end]

        # get its possible connections (if empty: start again from first root segment to go in other direction)
        connections = E₂(rg, e_current)
        if isempty(connections)
            if has_restarted
                fullgrown = true
            else
                has_restarted = true
            end
            continue #!
        end

        # filter on:
        # - remaining amount 
        # - edge2 != previous root segment
        # - edge2 has a different direction from edge1
        

        cidx = findfirst(
            is_valid_connection,
            connections
        )

        # pick a random edge2: 
            # add to root: 
                # if first run, append to end of root
                # if second run (going in other direction), append to start of root
            # remove root from remaining roots in that direction
            # remove connection
        # make new edge1 and repeat from 2
    end

end

function is_valid_connection(c, has_restarted, c_counts, es)
    if has_restarted
        (c_counts[c] > 0) &&
        !(es[2] ∈ c) &&
        follows_polarity[1] ? only(c[c != e_current]) : missing
    else
        c -> (c_counts[c] > 0) && !(es[end-1] ∈ c) && missing
    end
end
