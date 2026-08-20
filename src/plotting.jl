# PreGraph plotting
Plots.plot(bps::Vector{<:BranchPoint}; kwargs...) = scatter(x.(bps), y.(bps); kwargs...)
Plots.plot!(bps::Vector{<:BranchPoint}; kwargs...) = scatter!(x.(bps), y.(bps); kwargs...)

Plots.plot(bp::BranchPoint; kwargs...) = plot([bp]; kwargs...)
Plots.plot!(bp::BranchPoint; kwargs...) = plot!([bp]; kwargs...)

Plots.plot(g::PreGraph, s::Segment; kwargs...) = plot(
    [x(getbranchpoint(g, v)) for v in vertices(s)], [y(getbranchpoint(g, v)) for v in vertices(s)]; kwargs...
)
Plots.plot!(g::PreGraph, s::Segment; kwargs...) = plot!(
    [x(getbranchpoint(g, v)) for v in vertices(s)], [y(getbranchpoint(g, v)) for v in vertices(s)]; kwargs...
)

Plots.plot(g::PreGraph, ss::Vector{Segment{T, U}}; kwargs...) where {T, U} = plot(
    [x(getbranchpoint(g, v)) for s in ss for v in vertices(s)] |> x -> reshape(x, 2, :),
    [y(getbranchpoint(g, v)) for s in ss for v in vertices(s)] |> x -> reshape(x, 2, :);
    kwargs...
)

Plots.plot!(g::PreGraph, ss::Vector{Segment{T, U}}; kwargs...) where {T, U} = plot!(
    [x(getbranchpoint(g, v)) for s in ss for v in vertices(s)] |> x -> reshape(x, 2, :),
    [y(getbranchpoint(g, v)) for s in ss for v in vertices(s)] |> x -> reshape(x, 2, :);
    kwargs...
)

function Plots.plot(g::PreGraph; vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    p = plot(; kwargs...)
    plot!(g, segments(g); color = :black, edge_kwargs..., kwargs...)
    plot!(getbranchpoints(g); color = :grey, markersize = 2, vertex_kwargs..., kwargs...)
    return p
end

# RootGraph plotting

# vertices
Plots.plot(rvs::Vector{<:RootVertex}; kwargs...) = scatter(x.(rvs), y.(rvs); kwargs...)
Plots.plot!(rvs::Vector{<:RootVertex}; kwargs...) = scatter!(x.(rvs), y.(rvs); kwargs...)
Plots.plot(rv::RootVertex; kwargs...) = plot([rv]; kwargs...)
Plots.plot!(rv::RootVertex; kwargs...) = plot!([rv]; kwargs...)

# edges
Plots.plot(rg::RootGraph, res::Vector{<:RootEdge}; kwargs...) = plot(
    [x(getvertex(rg, v)) for re in res for v in vertices(re)] |> x -> reshape(x, 2, :),
    [y(getvertex(rg, v)) for re in res for v in vertices(re)] |> x -> reshape(x, 2, :);
    kwargs...
)
Plots.plot!(rg::RootGraph, res::Vector{<:RootEdge}; kwargs...) = plot!(
    [x(getvertex(rg, v)) for re in res for v in vertices(re)] |> x -> reshape(x, 2, :),
    [y(getvertex(rg, v)) for re in res for v in vertices(re)] |> x -> reshape(x, 2, :);
    kwargs...
)
Plots.plot(rg::RootGraph, re::RootEdge; kwargs...) = plot(rg, [re]; kwargs...)
Plots.plot!(rg::RootGraph, re::RootEdge; kwargs...) = plot!(rg, [re]; kwargs...)

# segments (one line per segment of the scan, using a representative edge)
Plots.plot(rg::RootGraph, segs::Vector{<:Vector{<:RootEdge}}; kwargs...) = plot(rg, first.(segs); kwargs...)
Plots.plot!(rg::RootGraph, segs::Vector{<:Vector{<:RootEdge}}; kwargs...) = plot!(rg, first.(segs); kwargs...)

# graph
function Plots.plot!(rg::RootGraph; annotate::Bool = false, augmented_alpha = 0.1, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    plot!(
        rg, segments(rg); color = :black,
        alpha = [is_augmented(seg) ? augmented_alpha : 1.0 for seg in segments(rg)] |> x -> reshape(x, 1, :),
        edge_kwargs..., kwargs...
    )

    annotate && merge!(vertex_kwargs, Dict(:annotation => [(x(rv), y(rv), text("$(id(rv))", 8, :right, :bottom)) for rv in V(rg)]))
    return plot!(V(rg); color = :grey, markersize = 2, vertex_kwargs..., kwargs...)
end

function Plots.plot(rg::RootGraph; aspect_ratio = :equal, legend = false, annotate::Bool = false, augmented_alpha = 0.1, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    p = plot(; aspect_ratio, legend, kwargs...)
    plot!(rg; annotate, augmented_alpha, vertex_kwargs, edge_kwargs, kwargs...)
    return p
end

# multiple graphs
function Plots.plot(rgs::Vector{<:RootGraph}; augmented_alpha = 0.1, aspect_ratio = :equal, legend = false, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    p = plot(; aspect_ratio, legend, kwargs...)
    for rg in rgs
        plot!(rg; augmented_alpha, vertex_kwargs, edge_kwargs, kwargs...)
    end
    return p
end

# with classification
# ## segments
function Plots.plot!(rg::RootGraph, segment_classification_dict::Dict{<:Vector{<:RootEdge}, <:Complex}; annotate::Bool = false, augmented_alpha = 0.1, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    classification_edge_kwargs = Dict(
        :color => [(real(segment_classification_dict[seg]) > 0) * RGB(1.0, 0, 0) + (imag(segment_classification_dict[seg]) > 0) * RGB(0, 0, 1.0) for seg in segments(rg)] |> x -> reshape(x, 1, :),
        :linewidth => [abs(segment_classification_dict[seg]) for seg in segments(rg)] |> x -> reshape(x, 1, :),
    )
    classification_vertex_kwargs = Dict(
        :color => :grey,
        :markersize => 2,
    )

    merge!(edge_kwargs, classification_edge_kwargs)
    merge!(vertex_kwargs, classification_vertex_kwargs)
    annotate && merge!(vertex_kwargs, Dict(:annotation => [(x(rv), y(rv), text("$(id(rv))", 8, :right, :bottom)) for rv in V(rg)]))

    p = plot!(rg; augmented_alpha, edge_kwargs, vertex_kwargs, kwargs...)

    return p
end

function Plots.plot(rg::RootGraph, segment_classification_dict::Dict{<:Vector{<:RootEdge}, <:Complex}; annotate::Bool = false, augmented_alpha = 0.1, aspect_ratio = :equal, legend = false, vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...)
    p = plot(; aspect_ratio, legend, kwargs...)
    plot!(rg, segment_classification_dict; annotate, augmented_alpha, vertex_kwargs, edge_kwargs, kwargs...)

    return p
end

function Plots.plot(
        rgs::Vector{RootGraph{T, U}}, segment_classification_dict::Dict{Int64, <:Dict{<:Vector{<:RootEdge}, <:Complex}};
        annotate::Bool = false, augmented_alpha = 0.1, aspect_ratio = :equal, legend = false,
        vertex_kwargs = Dict([]), edge_kwargs = Dict([]), kwargs...
    ) where {T, U}

    p = plot(; aspect_ratio, legend, kwargs...)
    for (i, rg) in enumerate(rgs)
        plot!(rg, segment_classification_dict[i]; annotate, augmented_alpha, vertex_kwargs, edge_kwargs, kwargs...)
    end

    return p
end

# ## edges (aggregated per segment)
function Plots.plot!(rg::RootGraph, edge_classification_dict::Dict{<:RootEdge, <:Complex}; kwargs...)
    segment_classification_dict = Dict([seg => sum([edge_classification_dict[e] for e in seg]) for seg in segments(rg)])
    return plot!(rg, segment_classification_dict; kwargs...)
end

function Plots.plot(rg::RootGraph, edge_classification_dict::Dict{<:RootEdge, <:Complex}; kwargs...)
    segment_classification_dict = Dict([seg => sum([edge_classification_dict[e] for e in seg]) for seg in segments(rg)])
    return plot(rg, segment_classification_dict; kwargs...)
end


# Root plotting
function Plots.plot(r::Root; kwargs...)
    linestyle = is_primary(r) ? :solid : :dot
    return plot(xs(r), ys(r); linestyle, color = :black, aspect_ratio = :equal, linewidth = 2, kwargs...)
end
function Plots.plot!(r::Root; kwargs...)
    linestyle = is_primary(r) ? :solid : :dot
    return plot!(xs(r), ys(r); linestyle, color = :black, aspect_ratio = :equal, linewidth = 2, kwargs...)
end

function Plots.plot(rs::Vector{<:Root}; kwargs...)
    plot()
    for (i, r) in enumerate(rs)
        linestyle = is_primary(r) ? :solid : :dot
        color = is_primary(r) ? HSV(0, 1, 0) : HSV(range(0, 360, length = length(rs))[i], 1, 0.75)
        plot!(xs(r), ys(r); linestyle, color, label = "$i", linewidth = 2, kwargs...)
    end
    return plot!(aspect_ratio = :equal; kwargs...)
end
function Plots.plot!(rs::Vector{<:Root}; kwargs...)
    for (i, r) in enumerate(rs)
        linestyle = is_primary(r) ? :solid : :dot
        color = is_primary(r) ? HSV(0, 1, 0) : HSV(range(0, 360, length = length(rs))[i], 1, 0.75)
        plot!(xs(r), ys(r); linestyle, color, label = "$i", linewidth = 2, kwargs...)
    end
    return
end

function Plots.plot(rss::Vector{<:Vector{<:Root}}; kwargs...)
    plot()
    for rs in rss
        plot!(rs; kwargs...)
    end
    return plot!(aspect_ratio = :equal; kwargs...)
end
function Plots.plot!(rss::Vector{<:Vector{<:Root}}; kwargs...)
    for rs in rss
        plot!(rs; kwargs...)
    end
    return
end

# custom plots
"""
    hypothesis_plot(rg::RootGraph)

Visualise the maximum allowed number of roots per segment of a graph.
"""
function hypothesis_plot(rg::RootGraph)
    nₕ(seg) = max(length(srcs(seg)), length(dsts(seg))) # maximum number of roots on a segment
    nₕs = nₕ.(segments(rg))
    ΔH = 360 / (maximum(nₕs) + 1)

    plot(
        rg, size = (600, 800), label = false, edge_kwargs =
            Dict(
            :linewidth => [width(seg) / 2 for seg in segments(rg)] |> x -> reshape(x, 1, :),
            :linecolor => [HSV(ΔH * n, 1, 0.75) for n in nₕs] |> x -> reshape(x, 1, :),
        ),
    )
    plot!(
        fill(missing, 1, length(unique(nₕs))), fill(missing, 1, length(unique(nₕs))),
        lw = 2, legend = true, legendfontsize = 12,
        label = reshape(sort(unique(nₕs)), 1, :),
        linecolor = reshape([HSV(ΔH * n, 1, 0.75) for n in sort(unique(nₕs))], 1, :)
    )
    return plot!(
        rg, filter(seg -> nₕ(seg) == 1, segments₀(rg)), aspect_ratio = :equal,
        label = false, color = HSV(ΔH, 1, 0.75), lw = 5
    )
end
