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
alpha(he::HyperEdge, augmented_alpha) = is_augmented(he) ? augmented_alpha : 1.0

# ### No classification
function graphplot(
        sg::SuperGraph; augmented_alpha = 0.05, size = (600, 400),
        vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )

    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    color = [
        # edge is 2 vertices + hidden NaN vertex
        fill(RGBAf(0, 0, 0, alpha(he, augmented_alpha)), 3)
        for he in Eₕ(sg)
    ] |> x -> reduce(vcat, x)

    lines!(ax, sg, Eₕ(sg); color, edge_kwargs...)
    scatter!(ax, Vₕ(sg); color = :grey, markersize = 8, vertex_kwargs...)

    return f
end

# ### With hyperedge classification
function graphplot(
        sg::SuperGraph, he_classification_dict::Dict{<:HyperEdge, <:Complex};
        augmented_alpha = 0.1, size = (600, 400), vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    )

    classification_edge_kwargs = Dict(
        :color => [
            fill(
                (real(he_classification_dict[he]) > 0) * RGBAf(1.0, 0, 0, alpha(he, augmented_alpha)) + 
                    (imag(he_classification_dict[he]) > 0) * RGBAf(0, 0, 1.0, alpha(he, augmented_alpha)),
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

# ### With num_hypotheses classification
"""
    hypothesis_plot(sg::SuperGraph)

Visualise the maximum allowed number of roots per segment of a graph.
"""
function hypothesis_plot(sg::SuperGraph; size = (600, 400), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    nₕs = [
        maximum([length(vertices(hv)) for hv in gethypervertex.([Vₕ(sg)], vertices(he))])
            for he in Eₕ₀(sg)
    ]
    ΔH = 360 / (maximum(nₕs) + 1)

    hue(i, ΔH) = ΔH * i
    edge_kwargs = Dict(
        :linewidth => [
            fill(nₕ == 1 ? 3 : width(he) / 2, 3)
            for (he, nₕ) in zip(Eₕ₀(sg), nₕs)
        ] |> x -> reduce(vcat, x),
        :color => [
            fill(HSVA(hue(nₕ, ΔH), 1, 0.75, alpha(he, 0.0)), 3)
            for (he, nₕ) in zip(Eₕ₀(sg), nₕs)
        ] |> x -> reduce(vcat, x),
    )

    lines!(ax, sg, Eₕ₀(sg); edge_kwargs...)
    label_lines = [lines!(ax, NaN, NaN, color = HSV(hue(i, ΔH), 1, 0.75)) for i in 1:maximum(nₕs)]
    Legend(
        f[1, 2],
        label_lines,
        string.(1:maximum(nₕs))
    )

    return f
end


# ### With annotation classification
function annotation_plot(sg::SuperGraph; size = (600, 400), fontsize = 6, edge_kwargs::Dict = Dict(), kwargs...)

    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    color(he) = ismissing(pred_primary(he)) ? HSV(0, 1, 0) : HSV(200, 1, pred_primary(he))

    for he in Eₕ₀(sg)
        lines!(ax, sg, he; color = color(he), edge_kwargs...)
    end

    annotation_coords = [(mean(xs(sg, he)), mean(ys(sg, he))) for he in Eₕ₀(sg)]
    annotation_texts = [string(segment_id(he)) for he in Eₕ₀(sg)]
    annotation!(ax, annotation_coords; text = annotation_texts, shrink = (0, 0), color = :red, fontsize)

    return f
end

function annotation_plot(
        sg::SuperGraph, annotation_dict::Dict{<:HyperEdge, Vector{<:Integer}}; 
        size = (600, 400), edge_kwargs::Dict = Dict(), kwargs...
    )

    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    color(he, ann_dict) = isempty(ann_dict[he]) ? RGBf(0.75, 0.75, 0.75) :
        (1 ∈ ann_dict[he] ? RGBf(1.0, 0.0, 0.0) : RGBf(0.0, 0.0, 1.0))

    # annotation_edge_kwargs = Dict(
    #     :color => [
    #         fill(color(he, annotation_dict), 3)
    #         for he in Eₕ₀(sg)
    #     ] |> x -> reduce(vcat, x),
    # )
    # edge_kwargs = merge(edge_kwargs, annotation_edge_kwargs)
    # lines!(ax, sg, Eₕ₀(sg); edge_kwargs...)

    for he in Eₕ₀(sg)
        lines!(ax, sg, he; color = color(he, annotation_dict), edge_kwargs...)
    end

    annotation_coords = [(mean(xs(sg, he)), mean(ys(sg, he))) for he in Eₕ₀(sg) if !isempty(annotation_dict[he])]
    annotation_texts = [string(annotation_dict[he])[2:end-1] for he in Eₕ₀(sg) if !isempty(annotation_dict[he])]
    annotation!(ax, annotation_coords, text = annotation_texts)

    return f
end

function correct_annotation!(annotation_dict::Dict{<:HyperEdge, Vector{<:Integer}}, sg::SuperGraph)
    f = annotation_plot(sg, annotation_dict)
    on(events(f).mouseposition) do event
        println(event)
    end

end

# # Root plotting
# ## Type recipe for single roots
Makie.convert_arguments(::Type{<:Lines}, r::Root) = (xs(r), ys(r))

# ## Root systems
function rootplot!(ax::Makie.Axis, rs::Vector{<:Root}; kwargs...)
    for (i, r) in enumerate(rs)
        linestyle = is_primary(r) ? :solid : :dot
        color = is_primary(r) ? HSV(0, 1, 0) : HSV(range(0, 360, length = length(rs))[i], 1, 0.75)
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
function rootplot(rss::Vector{<:Vector{<:Root}}; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    for rs in rss
        rootplot!(ax, rs; line_kwargs...)
    end

    return f
end