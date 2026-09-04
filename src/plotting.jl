# # SuperGraph plotting

# ## Type recipes for vertices and edges
# ### Vertices
Makie.convert_arguments(::Type{<:Scatter}, av::AbstractVertex) = (x(av), y(av))
Makie.convert_arguments(::Type{<:Scatter}, avs::Vector{<:AbstractVertex}) = (x.(avs), y.(avs))

# ### Edges
Makie.convert_arguments(::Type{<:Lines}, sg::SuperGraph, ae::AbstractEdge) = (xs(sg, ae), ys(sg, ae))
Makie.convert_arguments(::Type{<:Lines}, sg::SuperGraph, aes::Vector{<:AbstractEdge}) = (
    reduce(vcat, [[xs(sg, ae); NaN] for ae in aes]), reduce(vcat, [[ys(sg, ae); NaN] for ae in aes])
)


# ## SuperGraph #! TODO: use Makie's @recipe
alpha(he::HyperEdge, standard_alpha = 1.0, augmented_alpha = 0.05) = is_augmented(he) ? augmented_alpha : standard_alpha

# ### No classification
function graphplot!(
        ax::Makie.Axis, sg::SuperGraph; standard_alpha = 1.0, 
        augmented_alpha = 0.05, vertex_kwargs = Dict([]), edge_kwargs = Dict([])
    )
    color = [
        # edge is 2 vertices + hidden NaN vertex
        fill(RGBAf(0, 0, 0, alpha(he, standard_alpha, augmented_alpha)), 3)
        for he in Eₕ(sg)
    ] |> x -> reduce(vcat, x)

    lines!(ax, sg, Eₕ(sg); color, edge_kwargs...)
    scatter!(ax, Vₕ(sg); color = :grey, markersize = 8, vertex_kwargs...)
end

function graphplot(
        sg::SuperGraph; standard_alpha = 1.0, augmented_alpha = 0.05, 
        size = (600, 400), vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    graphplot!(ax, sg; standard_alpha, augmented_alpha, vertex_kwargs, edge_kwargs)

    return f
end

# ### With hyperedge classification
function graphplot(
        sg::SuperGraph, he_classification_dict::Dict{<:HyperEdge, <:Complex}; standard_alpha = 1.0,
        augmented_alpha = 0.1, size = (600, 400), vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )

    classification_edge_kwargs = Dict(
        :color => [
            fill(
                (real(he_classification_dict[he]) > 0) * RGBAf(1.0, 0, 0, alpha(he, standard_alpha, augmented_alpha)) + 
                    (imag(he_classification_dict[he]) > 0) * RGBAf(0, 0, 1.0, alpha(he, standard_alpha, augmented_alpha)),
                3
            )
            for he in Eₕ(sg)
        ] |> x -> reduce(vcat, x),
        :linewidth => [
            fill(abs(he_classification_dict[he]), 3) for he in Eₕ(sg)
        ] |> x -> reduce(vcat, x)
    )
    classification_vertex_kwargs = Dict(
        :color => :grey,
        :markersize => 8,
    )

    edge_kwargs = merge(edge_kwargs, classification_edge_kwargs)
    vertex_kwargs = merge(vertex_kwargs, classification_vertex_kwargs)

    return graphplot(sg; augmented_alpha, size, edge_kwargs, vertex_kwargs, kwargs...)
end

# ### With annotation classification
function add_annotation!(ax::Makie.Axis, sg::SuperGraph, hes::Vector{<:HyperEdge}; fontsize)
    annotation_coords = [(mean(xs(sg, he)), mean(ys(sg, he))) for he in hes]
    annotation_texts = [string(segment_id(he)) for he in hes]
    annotation!(ax, annotation_coords; text = annotation_texts, shrink = (0, 0), color = :red, fontsize)

    return nothing
end

function annotation_plot(sg::SuperGraph; size = (600, 400), fontsize = 6, edge_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    color(he) = ismissing(pred_primary(he)) ? HSV(0, 1, 0) : HSV(200, 1, pred_primary(he))

    for he in Eₕ₀(sg)
        lines!(ax, sg, he; color = color(he), edge_kwargs...)
    end

    add_annotation!(ax, sg, Eₕ₀(sg); fontsize)

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
function annotation_plot(sg::SuperGraph, root_system::Vector{<:Vector{<:Root}};
        size = (600, 400), fontsize = 6, edge_kwargs::Dict = Dict(), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    rootplot!(ax, root_system; edge_kwargs...)
    add_annotation!(ax, sg, Eₕ₀(sg); fontsize)

    return f
end

function annotation_plot(sg::SuperGraph, root_systems::Vector{<:Vector{<:Vector{<:Root}}};
        size = (600, 400), fontsize = 6, edge_kwargs::Dict = Dict(), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    graphplot!(
        ax, sg, standard_alpha = 0.5, augmented_alpha = 0.0, vertex_kwargs = Dict(:markersize => 3.0)
    )
    
    for root_system in root_systems
        rootplot!(ax, root_system; edge_kwargs...)
        add_annotation!(ax, sg, Eₕ₀(sg); fontsize)
    end

    return f
end

function annotation_plot(sgs::Vector{<:SuperGraph}, root_systems::Vector{<:Vector{<:Vector{<:Root}}};
        size = (600, 400), fontsize = 6, edge_kwargs::Dict = Dict(), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    
    for (sg, root_system) in zip(sgs, root_systems)
        graphplot!(
            ax, sg, standard_alpha = 0.5, augmented_alpha = 0.0, vertex_kwargs = Dict(:markersize => 3.0)
        )
        rootplot!(ax, root_system; edge_kwargs...)
        add_annotation!(ax, sg, Eₕ₀(sg); fontsize)
    end

    return f
end