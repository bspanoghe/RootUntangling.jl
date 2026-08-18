# PreGraph plotting
Plots.plot(vs::Vector{<:MetaVertex}; kwargs...) = scatter(x.(vs), y.(vs); kwargs...)
Plots.plot!(vs::Vector{<:MetaVertex}; kwargs...) = scatter!(x.(vs), y.(vs); kwargs...)

Plots.plot(v::MetaVertex; kwargs...) = plot([v]; kwargs...)
Plots.plot!(v::MetaVertex; kwargs...) = plot!([v]; kwargs...)

Plots.plot(g::PreGraph, s::Segment; kwargs...) = plot(
    [x(getmetavertex(g, v)) for v in vertices(s)], [y(getmetavertex(g, v)) for v in vertices(s)]; kwargs...
)
Plots.plot!(g::PreGraph, s::Segment; kwargs...) = plot!(
    [x(getmetavertex(g, v)) for v in vertices(s)], [y(getmetavertex(g, v)) for v in vertices(s)]; kwargs...
)

Plots.plot(g::PreGraph, ss::Vector{Segment{T, U}}; kwargs...) where {T, U} = plot(
    [x(getmetavertex(g, v)) for s in ss for v in vertices(s)] |> x -> reshape(x, 2, :),
    [y(getmetavertex(g, v)) for s in ss for v in vertices(s)] |> x -> reshape(x, 2, :);
    kwargs...
)

Plots.plot!(g::PreGraph, ss::Vector{Segment{T, U}}; kwargs...) where {T, U} = plot!(
    [x(getmetavertex(g, v)) for s in ss for v in vertices(s)] |> x -> reshape(x, 2, :),
    [y(getmetavertex(g, v)) for s in ss for v in vertices(s)] |> x -> reshape(x, 2, :);
    kwargs...
)

function Plots.plot(g::PreGraph; vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    p = plot(; kwargs...)
    plot!(g, segments(g); color = :black, edge_kwargs..., kwargs...)
    plot!(getmetavertices(g); color = :grey, markersize = 2, vertex_kwargs..., kwargs...)
    return p
end

# SuperGraph plotting

# vertices
Plots.plot(avs::Vector{<:AbstractVertex}; kwargs...) = scatter(x.(avs), y.(avs); kwargs...)
Plots.plot!(avs::Vector{<:AbstractVertex}; kwargs...) = scatter!(x.(avs), y.(avs); kwargs...)
Plots.plot(av::AbstractVertex; kwargs...) = plot([av]; kwargs...)
Plots.plot!(av::AbstractVertex; kwargs...) = plot!([av]; kwargs...)

# edges
Plots.plot(sg::SuperGraph, ses::Vector{<:SingularEdge}; kwargs...) = plot(
    [x(getsingularvertex(sg, v)) for se in ses for v in vertices(se)] |> x -> reshape(x, 2, :),
    [y(getsingularvertex(sg, v)) for se in ses for v in vertices(se)] |> x -> reshape(x, 2, :);
    kwargs...
)
Plots.plot!(sg::SuperGraph, ses::Vector{<:SingularEdge}; kwargs...) = plot!(
    [x(getsingularvertex(sg, v)) for se in ses for v in vertices(se)] |> x -> reshape(x, 2, :),
    [y(getsingularvertex(sg, v)) for se in ses for v in vertices(se)] |> x -> reshape(x, 2, :);
    kwargs...
)
Plots.plot(sg::SuperGraph, se::SingularEdge; kwargs...) = plot(sg, [se]; kwargs...)
Plots.plot!(sg::SuperGraph, se::SingularEdge; kwargs...) = plot!(sg, [se]; kwargs...)

Plots.plot(sg::SuperGraph, hes::Vector{<:HyperEdge}; kwargs...) = plot(
    [x(gethypervertex(sg, v)) for he in hes for v in vertices(he)] |> x -> reshape(x, 2, :),
    [y(gethypervertex(sg, v)) for he in hes for v in vertices(he)] |> x -> reshape(x, 2, :);
    kwargs...
)
Plots.plot!(sg::SuperGraph, hes::Vector{<:HyperEdge}; kwargs...) = plot!(
    [x(gethypervertex(sg, v)) for he in hes for v in vertices(he)] |> x -> reshape(x, 2, :),
    [y(gethypervertex(sg, v)) for he in hes for v in vertices(he)] |> x -> reshape(x, 2, :);
    kwargs...
)
Plots.plot(sg::SuperGraph, he::HyperEdge; kwargs...) = plot(sg, [he]; kwargs...)
Plots.plot!(sg::SuperGraph, he::HyperEdge; kwargs...) = plot!(sg, [he]; kwargs...)

# graph
function Plots.plot!(sg::SuperGraph; annotate::Bool = false, augmented_alpha = 0.1, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    plot!(
        sg, Eₕ(sg); color = :black, 
        alpha = [is_augmented(he) ? augmented_alpha : 1.0 for he in Eₕ(sg)] |> x -> reshape(x, 1, :),
        edge_kwargs..., kwargs...
    )

    annotate && merge!(vertex_kwargs, Dict(:annotation => [(x(hv), y(hv), text("$(id(hv))", 8, :right, :bottom)) for hv in Vₕ(sg)]))
    plot!(Vₕ(sg); color = :grey, markersize = 2, vertex_kwargs..., kwargs...)
end

function Plots.plot(sg::SuperGraph; aspect_ratio = :equal, legend = false, annotate::Bool = false, augmented_alpha = 0.1, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    p = plot(; aspect_ratio, legend, kwargs...)
    plot!(sg; annotate, augmented_alpha, vertex_kwargs, edge_kwargs, kwargs...)
    return p
end

# multiple graphs
function Plots.plot(sgs::Vector{<:SuperGraph}; augmented_alpha = 0.1, aspect_ratio = :equal, legend = false, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    p = plot(; aspect_ratio, legend, kwargs...)
    for sg in sgs
        plot!(sg; augmented_alpha, vertex_kwargs, edge_kwargs, kwargs...)
    end
    return p
end

# with classification
# ## hyperedges
function Plots.plot!(sg::SuperGraph, he_classification_dict::Dict{<:HyperEdge, <:Complex}; annotate::Bool = false, augmented_alpha = 0.1, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    classification_edge_kwargs = Dict(
        :color => [(real(he_classification_dict[he]) > 0) * RGB(1.0, 0, 0) + (imag(he_classification_dict[he]) > 0) * RGB(0, 0, 1.0) for he in Eₕ(sg)] |> x -> reshape(x, 1, :),
        :linewidth => [abs(he_classification_dict[he]) for he in Eₕ(sg)] |> x -> reshape(x, 1, :),
    )
    classification_vertex_kwargs = Dict(
        :color => :grey, 
        :markersize => 2,
    )

    merge!(edge_kwargs, classification_edge_kwargs)
    merge!(vertex_kwargs, classification_vertex_kwargs)
    annotate && merge!(vertex_kwargs, Dict(:annotation => [(x(hv), y(hv), text("$(id(hv))", 8, :right, :bottom)) for hv in Vₕ(sg)]))

    p = plot!(sg; augmented_alpha, edge_kwargs, vertex_kwargs, kwargs...)

    return p
end

function Plots.plot(sg::SuperGraph, he_classification_dict::Dict{<:HyperEdge, <:Complex}; annotate::Bool = false, augmented_alpha = 0.1, aspect_ratio = :equal, legend = false, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    p = plot(; aspect_ratio, legend, kwargs...)
    plot!(sg, he_classification_dict; annotate, augmented_alpha, vertex_kwargs, edge_kwargs, kwargs...)

    return p
end

function Plots.plot(sgs::Vector{SuperGraph{T, U}}, he_classification_dict::Dict{Int64, Dict{HyperEdge{T, U}, Complex{Int64}}};
        annotate::Bool = false, augmented_alpha = 0.1, aspect_ratio = :equal, legend = false,
        vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    ) where {T, U}
        
    p = plot(; aspect_ratio, legend, kwargs...)
    for (i, sg) in enumerate(sgs)
        plot!(sg, he_classification_dict[i]; annotate, augmented_alpha, vertex_kwargs, edge_kwargs, kwargs...)
    end

    return p
end

# ## singular edges
function Plots.plot!(sg::SuperGraph, se_classification_dict::Dict{<:SingularEdge, <:Complex}; annotate::Bool = false, augmented_alpha = 0.1, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    classification_edge_kwargs = Dict(
        :color => [(real(se_classification_dict[he]) > 0) * RGB(1.0, 0, 0) + (imag(se_classification_dict[he]) > 0) * RGB(0, 0, 1.0) for he in Eₕ(sg)] |> x -> reshape(x, 1, :),
        :linewidth => [abs(se_classification_dict[he]) for he in Eₕ(sg)] |> x -> reshape(x, 1, :),
    )
    classification_vertex_kwargs = Dict(
        :color => :grey, 
        :markersize => 2,
    )

    merge!(edge_kwargs, classification_edge_kwargs)
    merge!(vertex_kwargs, classification_vertex_kwargs)
    annotate && merge!(vertex_kwargs, Dict(:annotation => [(x(hv), y(hv), text("$(id(hv))", 8, :right, :bottom)) for hv in Vₕ(sg)]))

    p = plot!(sg; augmented_alpha, edge_kwargs, vertex_kwargs, kwargs...)

    return p
end

function Plots.plot(sg::SuperGraph, se_classification_dict::Dict{<:SingularEdge, <:Complex}; annotate::Bool = false, augmented_alpha = 0.1, aspect_ratio = :equal, legend = false, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    p = plot(; aspect_ratio, legend, kwargs...)
    plot!(sg, se_classification_dict; annotate, augmented_alpha, vertex_kwargs, edge_kwargs, kwargs...)

    return p
end


# Root plotting
function Plots.plot(r::Root; kwargs...)
    linestyle = is_primary(r) ? :solid : :dot
    plot(xs(r), ys(r); linestyle, color = :black, aspect_ratio = :equal, linewidth = 2, kwargs...)
end
function Plots.plot!(r::Root; kwargs...)
    linestyle = is_primary(r) ? :solid : :dot
    plot!(xs(r), ys(r); linestyle, color = :black, aspect_ratio = :equal, linewidth = 2, kwargs...)
end

function Plots.plot(rs::Vector{<:Root}; kwargs...)
    plot()
    for (i, r) in enumerate(rs)
        linestyle = is_primary(r) ? :solid : :dot
        color = is_primary(r) ? HSV(0, 1, 0) : HSV(range(0, 360, length = length(rs))[i], 1, 0.75)
        plot!(xs(r), ys(r); linestyle, color, label = "$i", linewidth = 2, kwargs...)
    end
    plot!(aspect_ratio = :equal; kwargs...)
end
function Plots.plot!(rs::Vector{<:Root}; kwargs...)
    for (i, r) in enumerate(rs)
        linestyle = is_primary(r) ? :solid : :dot
        color = is_primary(r) ? HSV(0, 1, 0) : HSV(range(0, 360, length = length(rs))[i], 1, 0.75)
        plot!(xs(r), ys(r); linestyle, color, label = "$i", linewidth = 2, kwargs...)
    end
end

function Plots.plot(rss::Vector{<:Vector{<:Root}}; kwargs...)
    plot()
    for rs in rss
        plot!(rs; kwargs...)
    end
    plot!(aspect_ratio = :equal; kwargs...)
end
function Plots.plot!(rss::Vector{<:Vector{<:Root}}; kwargs...)
    for rs in rss
        plot!(rs; kwargs...)
    end
end

# custom plots
"""
    hypothesis_plot(sg::SuperGraph)

Visualise the maximum allowed number of roots per segment of a graph.
"""
function hypothesis_plot(sg::SuperGraph)
    nₕs = [
        maximum([length(vertices(hv)) for hv in gethypervertex.([Vₕ(sg)], vertices(he))])
        for he in Eₕ(sg)
    ]
    nₕ₀s = [
        maximum([length(vertices(hv)) for hv in gethypervertex.([Vₕ(sg)], vertices(he))])
        for he in Eₕ₀(sg)
    ]
    ΔH = 360 / (maximum(nₕs) + 1)

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
    plot!(sg, HyperEdge[he for (nₕ₀, he) in zip(nₕ₀s, Eₕ₀(sg)) if nₕ₀ == 1], aspect_ratio = :equal,
        label = false, color = HSV(ΔH, 1, 0.75), lw = 5
    )
end