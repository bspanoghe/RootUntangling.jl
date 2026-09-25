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
        augmented_alpha = 0.1, size = (600, 400), vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )

    rd = get_result_dict(rg, model)

    classification_edge_kwargs = Dict(
        :color => [
            fill(
                rd[re][:ep] * RGBAf(1.0, 0, 0, alpha(re, standard_alpha, augmented_alpha)) + 
                    rd[re][:el] * RGBAf(0, 0, 1.0, alpha(re, standard_alpha, augmented_alpha)),
                3
            )
            for re in E(rg)
        ] |> x -> reduce(vcat, x),
        :linewidth => [
            fill(rd[re][:en], 3) for re in E(rg)
        ] |> x -> reduce(vcat, x)
    )
    classification_vertex_kwargs = Dict(:color => :grey)

    edge_kwargs = merge(edge_kwargs, classification_edge_kwargs)
    vertex_kwargs = merge(vertex_kwargs, classification_vertex_kwargs)

    return graphplot(rg; augmented_alpha, size, edge_kwargs, vertex_kwargs, kwargs...)
end

# ### With annotation classification
function add_annotation!(ax::Makie.Axis, rg::RootGraph, res::Vector{<:RootEdge}; fontsize)
    annotation_coords = [(mean(xs(rg, re)), mean(ys(rg, re))) for re in res]
    annotation_texts = [string(segment_id(re)) for re in res]
    annotation!(ax, annotation_coords; text = annotation_texts, shrink = (0, 0), color = :red, fontsize)

    return nothing
end

function annotation_plot(rg::RootGraph; size = (600, 400), fontsize = 6, edge_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    color(re) = ismissing(pred_primary(re)) ? HSV(0, 1, 0) : HSV(200, 1, pred_primary(re))

    for re in Eₕ₀(rg)
        lines!(ax, rg, re; color = color(re), edge_kwargs...)
    end

    add_annotation!(ax, rg, Eₕ₀(rg); fontsize)

    return f
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

# ### Multiple (entangled) root systems
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

# ### annotated rootplot
function annotation_plot(rg::RootGraph, root_system::Vector{<:Vector{<:Root}};
        size = (600, 400), fontsize = 6, edge_kwargs::Dict = Dict(), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    rootplot!(ax, root_system; edge_kwargs...)
    add_annotation!(ax, rg, Eₕ₀(rg); fontsize)

    return f
end

function annotation_plot(rg::RootGraph, root_systems::Vector{<:Vector{<:Vector{<:Root}}};
        size = (600, 400), fontsize = 6, edge_kwargs::Dict = Dict(), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    graphplot!(
        ax, rg, standard_alpha = 0.5, augmented_alpha = 0.0, vertex_kwargs = Dict(:markersize => 3.0)
    )
    
    for root_system in root_systems
        rootplot!(ax, root_system; edge_kwargs...)
        add_annotation!(ax, rg, Eₕ₀(rg); fontsize)
    end

    return f
end

function annotation_plot(rgs::Vector{<:RootGraph}, root_systems::Vector{<:Vector{<:Vector{<:Root}}};
        size = (600, 400), fontsize = 6, edge_kwargs::Dict = Dict(), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    
    for (rg, root_system) in zip(rgs, root_systems)
        graphplot!(
            ax, rg, standard_alpha = 0.5, augmented_alpha = 0.0, vertex_kwargs = Dict(:markersize => 3.0)
        )
        rootplot!(ax, root_system; edge_kwargs...)
        add_annotation!(ax, rg, Eₕ₀(rg); fontsize)
    end

    return f
end