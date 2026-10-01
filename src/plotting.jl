# # RootGraph plotting

# ## Type recipes for vertices and edges
# ### Vertices
Makie.convert_arguments(::Type{<:Scatter}, rv::RootVertex) = (x(rv), y(rv))
Makie.convert_arguments(::Type{<:Scatter}, rvs::Vector{<:RootVertex}) = (x.(rvs), y.(rvs))

# ### Edges
Makie.convert_arguments(::Type{<:Lines}, rg::RootGraph, re::RootEdge) = (xs(rg, re), ys(rg, re))
Makie.convert_arguments(::Type{<:Lines}, rg::RootGraph, res::Vector{<:RootEdge}) = (
    reduce(vcat, [[xs(rg, re); NaN] for re in res]), reduce(vcat, [[ys(rg, re); NaN] for re in res])
)


# ## RootGraph #! TODO: use Makie's @recipe
alpha(re::RootEdge, standard_alpha = 1.0, augmented_alpha = 0.05) = is_augmented(re) ? augmented_alpha : standard_alpha

# ### No classification
function graphplot!(
        ax::Makie.Axis, rg::RootGraph; standard_alpha = 1.0, 
        augmented_alpha = 0.05, vertex_kwargs = Dict([]), edge_kwargs = Dict([])
    )
    color = [
        # edge is 2 vertices + hidden NaN vertex
        fill(RGBAf(0, 0, 0, alpha(re, standard_alpha, augmented_alpha)), 3)
        for re in E(rg)
    ] |> x -> reduce(vcat, x)

    lines!(ax, rg, E(rg); color, edge_kwargs...)
    scatter!(ax, V(rg); color = :grey, markersize = 5, vertex_kwargs...)
end

function graphplot(
        rg::RootGraph; standard_alpha = 1.0, augmented_alpha = 0.05, 
        size = (600, 600), vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    graphplot!(ax, rg; standard_alpha, augmented_alpha, vertex_kwargs, edge_kwargs)

    return f
end

# ### With rootedge classification
function graphplot(
        rg::RootGraph, model::JuMP.Model; standard_alpha = 1.0,
        augmented_alpha = 0.1, size = (600, 600), vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )

    eps = round.(Int64, value.(model[:ep₊])) + round.(Int64, value.(model[:ep₋]))
    els = round.(Int64, value.(model[:el₊])) + round.(Int64, value.(model[:el₋]))

    classification_edge_kwargs = Dict(
        :color => [
            fill(
                (ep > 0) * RGBAf(1.0, 0, 0, alpha(re, standard_alpha, augmented_alpha)) + 
                    (el > 0) * RGBAf(0, 0, 1.0, alpha(re, standard_alpha, augmented_alpha)),
                3
            )
            for (re, ep, el) in zip(E(rg), eps, els)
        ] |> x -> reduce(vcat, x),
        :linewidth => [
            fill(ep + el, 3) for (ep, el) in zip(eps, els)
        ] |> x -> reduce(vcat, x)
    )
    classification_vertex_kwargs = Dict(:color => :grey)

    edge_kwargs = merge(edge_kwargs, classification_edge_kwargs)
    vertex_kwargs = merge(vertex_kwargs, classification_vertex_kwargs)

    return graphplot(rg; augmented_alpha, size, edge_kwargs, vertex_kwargs, kwargs...)
end

# ### With annotation classification
function annotation_plot(rg::RootGraph; size = (500, 1000), fontsize = 8, e_color = :red,
        edge_kwargs::Dict = Dict(), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    color(re) = ismissing(pred_primary(re)) ? HSV(0, 1, 0) : HSV(200, 1, pred_primary(re))

    for re in E₀(rg)
        lines!(ax, rg, re; color = color(re), edge_kwargs...)
    end

    annotation_coords = [(mean(xs(rg, re)), mean(ys(rg, re))) for re in E₀(rg)]
    annotation_texts = [string(segment_id(re)) for re in E₀(rg)]
    annotation!(ax, annotation_coords; text = annotation_texts, shrink = (0, 0), color = e_color, fontsize)

    return f
end

function annotation_plot(rg::RootGraph, model::JuMP.Model;
        fontsize = 8, e_color = :red, v_color = :green,
        standard_alpha = 1.0, augmented_alpha = 0.1, size = (500, 1000), 
        vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )

    f_g = graphplot(rg, model; standard_alpha, augmented_alpha, size, vertex_kwargs, edge_kwargs, kwargs...)

    e₊s = round.(Int64, value.(model[:ep₊])) + round.(Int64, value.(model[:el₊]))
    e₋s = round.(Int64, value.(model[:ep₋])) + round.(Int64, value.(model[:el₋]))

    segment_coords = [(mean(xs(rg, re)), mean(ys(rg, re))) for re in E₀(rg)]
    segment_texts = [
        "($e₊ / $e₋)"
        for (re, e₊, e₋) in zip(E(rg), e₊s, e₋s)
        if !is_augmented(re)
    ]
    annotation!(f_g.content[1], segment_coords; text = segment_texts, color = e_color, shrink = (0, 0), fontsize)

    vertex_coords = [(x(rv), y(rv)) for rv in V₀(rg)]
    vertex_texts = string.(id.(V₀(rg)))
    annotation!(f_g.content[1], vertex_coords; text = vertex_texts, color = v_color, shrink = (0, 0), fontsize)
    return f_g
end

# # Root plotting
# ## Type recipe for single roots
Makie.convert_arguments(::Type{<:Lines}, rg::RootGraph, r::Root) = (xs(rg, r), ys(rg, r))

# ## Root systems
function rootplot!(ax::Makie.Axis, rg::RootGraph, rs::RootSystem; kwargs...)
    r = primary(rs)
    lines!(ax, rg, r; color = HSV(0, 1, 0), linestyle = :solid, label = "0", linewidth = 2, kwargs...)

    for (i, r) in enumerate(laterals(rs))
        color = HSV(range(0, 360, length = length(rs))[i], 1, 0.75)
        lines!(ax, rg, r; color, linestyle = :dot, label = "$i", linewidth = 2, kwargs...)
    end
end

function rootplot(rg::RootGraph, rs::RootSystem; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, rg, rs; line_kwargs...)

    return f
end

function rootplot!(ax::Makie.Axis, rg::RootGraph, rss::Vector{<:RootSystem}; line_kwargs::Dict = Dict(), kwargs...)
    for rs in rss
        rootplot!(ax, rg, rs; line_kwargs...)
    end
end

function rootplot(rg::RootGraph, rss::Vector{<:RootSystem}; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, rg, rss; line_kwargs...)

    return f
end