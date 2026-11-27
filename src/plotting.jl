# # RootGraph plotting

# ## Type recipes for vertices and edges
# ### Vertices
Makie.convert_arguments(::Type{<:Scatter}, rv::RootVertex) = (x(rv), y(rv))
Makie.convert_arguments(::Type{<:Scatter}, rvs::Vector{<:RootVertex}) = (x.(rvs), y.(rvs))

# ### Edges
Makie.convert_arguments(::Type{<:Lines}, rg::RootGraph, ae::AbstractEdge) = (xs(rg, ae), ys(rg, ae))
Makie.convert_arguments(::Type{<:Lines}, rg::RootGraph, aes::Vector{<:AbstractEdge}) = (
    reduce(vcat, [[xs(rg, ae); NaN] for ae in aes]), reduce(vcat, [[ys(rg, ae); NaN] for ae in aes])
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
        v_color = :green, arrow_scaling::Real = 20, use_plain_text::Bool = false, edge_kwargs::Dict = Dict(), kwargs...
    )
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)

    color(re) = ismissing(pred_primary(re)) ? HSV(0, 1, 0) : HSV(0, 1, pred_primary(re))

    for re in E₀(rg)
        lines!(ax, rg, re; color = color(re), edge_kwargs...)
    end

    arrows2d!(ax,
        first.(xs.([rg], E₀(rg))), first.(ys.([rg], E₀(rg))), 
        cos.(src_angle.(E₀(rg))), sin.(src_angle.(E₀(rg))),
        lengthscale = arrow_scaling, color = :orange, shaftwidth = 2, tipwidth = 5, tiplength = 5
    )
    arrows2d!(ax,
        last.(xs.([rg], E₀(rg))), last.(ys.([rg], E₀(rg))),
        cos.(dst_angle.(E₀(rg))), sin.(dst_angle.(E₀(rg))),
        lengthscale = arrow_scaling, color = :cyan, shaftwidth = 2, tipwidth = 5, tiplength = 5
    )

    annotation_coords = [(mean(xs(rg, re)), mean(ys(rg, re))) for re in E₀(rg)]
    annotation_texts = [string(segment_id(re)) for re in E₀(rg)]
    
    vertex_coords = [(x(rv), y(rv)) for rv in V₀(rg)]
    vertex_texts = string.(id.(V₀(rg)))
    if use_plain_text
        text!(ax, annotation_coords; text = annotation_texts, color = e_color, fontsize)
        text!(ax, vertex_coords; text = vertex_texts, color = v_color, fontsize)
    else
        annotation!(ax, annotation_coords; text = annotation_texts, shrink = (0, 0), color = e_color, fontsize)
        annotation!(ax, vertex_coords; text = vertex_texts, color = v_color, shrink = (0, 0), fontsize)
    end

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
Makie.convert_arguments(::Type{<:Lines}, r::Root) = (xs(r), ys(r))
Makie.convert_arguments(::Type{<:Scatter}, vr::ValRoot; σ = 0.5) = (xs(vr) .+ σ*randn(length(xs(vr))), ys(vr) .+ σ*randn(length(xs(vr))))

# ## Root systems
function rootplot!(ax::Makie.Axis, rs::RootSystem; kwargs...)
    r = primary(rs)
    lines!(ax, r; color = HSV(0, 1, 0), linestyle = :solid, label = "0", linewidth = 2, kwargs...)

    for (i, r) in enumerate(laterals(rs))
        color = HSV(range(0, 360, length = length(rs))[i], 1, 0.75)
        lines!(ax, r; color, linestyle = :dot, label = "$i", linewidth = 2, kwargs...)
    end
end

function rootplot(rs::RootSystem; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, rs; line_kwargs...)

    return f
end

function rootplot!(ax::Makie.Axis, rss::Vector{<:RootSystem}; line_kwargs::Dict = Dict(), kwargs...)
    for rs in rss
        rootplot!(ax, rs; line_kwargs...)
    end
end

function rootplot(rss::Vector{<:RootSystem}; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, rss; line_kwargs...)

    return f
end

function rootplot!(ax::Makie.Axis, vrs::ValRootSystem; label_prefix = "", markersize = 3, alpha = 0.75, kwargs...)
    r = primary(vrs)
    scatter!(ax, r; color = HSV(0, 1, 0), label = label_prefix * "0", markerspace = :data, markersize, alpha, kwargs...)

    for (i, r) in enumerate(laterals(vrs))
        color = HSV(range(0, 360, length = length(vrs))[i], 1, 0.75)
        scatter!(ax, r; color, label = label_prefix * "$i", markersize, alpha, markerspace = :data, kwargs...)
    end
end

function rootplot!(ax::Makie.Axis, vrss::Vector{<:ValRootSystem}; kwargs...)
    for (i, vrs) in enumerate(vrss)
        rootplot!(ax, vrs; label_prefix = "$i", kwargs...)
    end
end

function rootplot(vrs::ValRootSystem; size = (600, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, vrs; line_kwargs...)

    return f
end

function rootplot(vrss::Vector{<:ValRootSystem}; size = (1200, 400), line_kwargs::Dict = Dict(), kwargs...)
    f = Figure(; size)
    ax = Axis(f[1, 1]; aspect = DataAspect(), kwargs...)
    rootplot!(ax, vrss; line_kwargs...)

    return f
end

function rootplot(vrs::ValRootSystem, rs::RootSystem; size = (1200, 600), kwargs...)
    f = Figure(; size)
    rootplot!(Axis(f[1, 1], aspect = DataAspect()), vrs)
    f[1, 2] = Legend(f, f.content[end])
    rootplot!(Axis(f[1, 3], aspect = DataAspect()), rs)
    f[1, 4] = Legend(f, f.content[end])
    return f
end

# # more validation things
function plot_mistakes(rg::RootGraph, rss::Vector{<:RootSystem}, vrss::Vector{<:ValRootSystem}, rs_mapping, root_mappings)
    rss, vrss = (deepcopy(rss), deepcopy(vrss))
    for (i, (rs, vrs)) in enumerate(zip(rss, vrss[rs_mapping]))
        add_dummy_roots!(rs, vrs)
        laterals(vrs) .= laterals(vrs)[root_mappings[i]];
    end

    incorrect_segments_ids = [
        [
            full_setdiff(standard_segments(laterals(rs)[i]), standard_segments(laterals(vrs)[i]))
            for i in eachindex(laterals(rs))
        ] |> x -> reduce(vcat, x)
        for (rs, vrs) in zip(rss, vrss[rs_mapping])
    ] |> x -> reduce(vcat, x)

    f = graphplot(rg)
    for re in filter(re -> re.og_segment_id ∈ incorrect_segments_ids, E₀(rg))
        lines!(f[1, 1], rg, re, color = :red, linewidth = 10, alpha = 0.3)
    end
    return f
end