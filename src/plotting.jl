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
    scatter!(ax, V(rg); color = :grey, markersize = 8, vertex_kwargs...)
end

function graphplot(
        rg::RootGraph; standard_alpha = 1.0, augmented_alpha = 0.05, 
        size = (600, 400), vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    graphplot!(ax, rg; standard_alpha, augmented_alpha, vertex_kwargs, edge_kwargs)

    return f
end

# ### With hyperedge classification
function graphplot(
        rg::RootGraph, re_classification_dict::Dict{<:RootEdge, <:Complex}; standard_alpha = 1.0,
        augmented_alpha = 0.1, size = (600, 400), vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )

    classification_edge_kwargs = Dict(
        :color => [
            fill(
                (real(re_classification_dict[re]) > 0) * RGBAf(1.0, 0, 0, alpha(re, standard_alpha, augmented_alpha)) + 
                    (imag(re_classification_dict[re]) > 0) * RGBAf(0, 0, 1.0, alpha(re, standard_alpha, augmented_alpha)),
                3
            )
            for re in Eₕ(rg)
        ] |> x -> reduce(vcat, x),
        :linewidth => [
            fill(abs(re_classification_dict[re]), 3) for re in Eₕ(rg)
        ] |> x -> reduce(vcat, x)
    )
    classification_vertex_kwargs = Dict(
        :color => :grey,
        :markersize => 8,
    )

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
Makie.convert_arguments(::Type{<:Lines}, r::Root) = (xs(r), ys(r))

# ## Root systems
function rootplot!(ax::Makie.Axis, rs::Vector{<:Root}; kwargs...)
    for (i, r) in enumerate(rs)
        linestyle = is_primary(r) ? :solid : :dot
        color = is_primary(r) ? HSV(0, 1, 0) : HSV(range(0, 360, length = length(rs)+1)[i], 1, 0.75)
        lines!(ax, r; linestyle, color, label = "$i", linewidth = 2, kwargs...)
    end
end

function rootplot(rs::Vector{<:Root}; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, rs; line_kwargs...)

    return f
end

# ### Multiple (entangled) root systems
function rootplot!(ax::Makie.Axis, rss::Vector{<:Vector{<:Root}}; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    for rs in rss
        rootplot!(ax, rs; line_kwargs...)
    end
end

function rootplot(rss::Vector{<:Vector{<:Root}}; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, rss; line_kwargs...)

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