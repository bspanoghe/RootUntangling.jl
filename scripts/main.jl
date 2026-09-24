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

directory = validation_dir
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
        rg; optimizer = Gurobi.Optimizer, time_limit = 60,
        num_roots = 1
    )

    annotate_that_thang = true
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

c_counts = Dict(E₂(rg) .=> round.(Int64, value.(model[:cp])))
e_counts = Dict(E(rg) .=> round.(Int64, value.(model[:ep])))
e₊_counts = Dict(E(rg) .=> round.(Int64, value.(model[:ep₊])))
e₋_counts = Dict(E(rg) .=> round.(Int64, value.(model[:ep]) - value.(model[:ep₊])))

primary_fragments = fragment(rg, true, c_counts, e_counts, e₊_counts, e₋_counts)
unambiguous_stitch_fragments!(primary_fragments)

c_counts = Dict(E₂(rg) .=> round.(Int64, value.(model[:cl])))
e_counts = Dict(E(rg) .=> round.(Int64, value.(model[:el])))
e₊_counts = Dict(E(rg) .=> round.(Int64, value.(model[:el₊])))
e₋_counts = Dict(E(rg) .=> round.(Int64, value.(model[:el]) - value.(model[:el₊])))

lateral_fragments = fragment(rg, true, c_counts, e_counts, e₊_counts, e₋_counts)
unambiguous_stitch_fragments!(lateral_fragments)






# assign laterals that split to their primary root

# assign laterals that appeared to their most probable primary root




# get primary root(s)

import RootUntangling: RootEdge, RootFragment, DirectedRootFragment, UndirectedRootFragment, edge_vertices

function fragment(rg::RootGraph{T, U}, is_primary::Bool,
        c_counts::Dict, e_counts::Dict, e₊_counts::Dict, e₋_counts::Dict
    )

    fragments = RootFragment[]

    for c in E₂(rg)
        v_shared = shared_vertex(c)

        if c_counts[c] == e_counts[c[1]] # first edge fully explains connection

            e = c[1]
            for i in 1:e₊_counts[e] # root follows direction of edge
                if dst(e) == v_shared # does edge go toward shared vertex?
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = false)
                else
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = true)
                end                    
            end

            for i in 1:e₋_counts[e] # root goes against direction of edge
                if dst(e) != v_shared # flip direction
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = false)
                else
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = true)
                end        
            end

        elseif c_counts[c] == e_counts[c[2]] # second edge fully explains connection

            e = c[2]
            for i in 1:e₊_counts[e] # root follows direction of edge
                if src(e) == v_shared # does edge go away from shared vertex?
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = false)
                else
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = true)
                end
            end

            for i in 1:e₋_counts[e] # root goes against direction of edge
                if src(e) != v_shared # flip direction
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = false)
                else
                    add_directed_fragment!(fragments, is_primary, c, v_shared, switch_edges = true)
                end
            end

        else # neither edge fully explains connection

            for i in 1:c_counts[c]
                push!(fragments, UndirectedRootFragment(is_primary, c))
            end

        end
    end

    return fragments
end

shared_vertex(c::Vector{<:RootEdge}) = vertices(c[1])[findfirst(v -> in(v, vertices(c[2])), vertices(c[1]))]
function add_directed_fragment!(fragments::Vector{<:RootFragment}, is_primary::Bool, c::Vector{<:RootEdge}, v_shared; switch_edges::Bool)
    if !switch_edges
        push!(fragments, DirectedRootFragment(
            is_primary, 
            [
                RootArc(c[1], keep_order = (dst(c[1]) == v_shared)),
                RootArc(c[2], keep_order = (src(c[2]) == v_shared))
            ])
        )
    else
        push!(fragments, DirectedRootFragment(is_primary,
            [
                RootArc(c[2], keep_order = (dst(c[2]) == v_shared)),
                RootArc(c[1], keep_order = (src(c[1]) == v_shared))
            ])
        )
    end

    return nothing
end



function unambiguous_stitch_fragments!(fragments::Vector{<:RootFragment})

    growing = true
    # keep going until all root fragments cant grow anymore
    while growing
        growing = false
        for (i, fragment) in enumerate(fragments)
            
            f_end_idxs = findall(f -> connects_to_end(fragment, f), fragments)
            if are_options_unambiguous(fragments[f_end_idxs])
                stitch_to_end!(fragment, fragments[f_end_idxs[1]])
                deleteat!(fragments, f_end_idxs[1])
                
                growing = true
            end

            f_start_idxs = findall(f -> connects_to_start(fragment, f), fragments)
            if are_options_unambiguous(fragments[f_start_idxs])
                stitch_to_start!(fragment, fragments[f_start_idxs[1]]) # connect one of them
                deleteat!(fragments, f_start_idxs[1])

                growing = true
            end

        end
    end

    return nothing
end

are_options_unambiguous(xs) = (
    (length(xs) == 1) || # there is only one option
        (length(xs) > 1 && allequal(xs)) # all options are equal
)

connects_to_end(f1::DirectedRootFragment, f2::DirectedRootFragment) = edge_vertices(f1)[end] == edge_vertices(f2)[1]
connects_to_start(f1::DirectedRootFragment, f2::DirectedRootFragment) = edge_vertices(f1)[1] == edge_vertices(f2)[end]

stitch_to_end!(f1::DirectedRootFragment, f2::DirectedRootFragment) = append!(f1.edge_vertices, f2.edge_vertices[2:end])
stitch_to_start!(f1::DirectedRootFragment, f2::DirectedRootFragment) = prepend!(f1.edge_vertices, f2.edge_vertices[1:end-1])

# https://www.youtube.com/watch?v=TVI7S6zKiBM
function ambiguous_stitch_fragments!(fragments::Vector{<:RootFragment})

end








function get_possible_directions(c::Vector{<:RootEdge}, e₊_counts::Dict, e₋_counts::Dict)
    opposite_direction = allequal(src.(c)) || allequal(dst.(c)) # both edges point towards or away from shared vertex
    return Dict(
        [true, true] => !opposite_direction && e₊_counts[c[1]] > 0 && e₊_counts[c[2]] > 0,
        [false, false] => !opposite_direction && e₋_counts[c[1]] > 0 && e₋_counts[c[2]] > 0,
        [true, false] => opposite_direction && e₊_counts[c[1]] > 0 && e₋_counts[c[2]] > 0,
        [false, true] => opposite_direction && e₋_counts[c[1]] > 0 && e₊_counts[c[2]] > 0,
    )
end

function grow_rootfragment!(c_counts::Dict, e₊_counts::Dict, e₋_counts::Dict, rg::RootGraph, root_type::Symbol)
    root_type ∈ [:primary, :lateral] || error("Root type should be primary or lateral")
    
    # choose directions of edges
    directions = findfirst(get_possible_directions(c0, e₊_counts, e₋_counts))
    isnothing(directions) && error("Oh what the heck")
    RootFragment(
        root_type == :primary,
        c0,
        vertices(c0[1])[findfirst(v -> !in(v, vertices(c0[2])), vertices(c0[1]))],
        vertices(c0[2])[findfirst(v -> !in(v, vertices(c0[1])), vertices(c0[2]))],
        directions...
    )
    
    # instantiate root fragment

    
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
