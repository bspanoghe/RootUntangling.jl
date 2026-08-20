# # classification of graph vertices/edges/segments
# the classification is a complex number: the real part counts primary roots, the imaginary part lateral roots
# one graph

function get_edge_classification_dict(rg::RootGraph, model::JuMP.Model)
    ea = value.(model[:ea])
    ep = value.(model[:ep])

    return Dict([E(rg)[i] => round(Bool, ep[i]) + round(Bool, ea[i] - ep[i]) * im for i in eachindex(E(rg))])
end

function get_segment_classification_dict(rg::RootGraph, model::JuMP.Model)
    edge_classification_dict = get_edge_classification_dict(rg, model)

    return Dict([seg => sum([edge_classification_dict[e] for e in seg]) for seg in segments(rg)])
end

function get_vertex_classification_dict(rg::RootGraph, model::JuMP.Model)
    va = value.(model[:va])
    vp = value.(model[:vp])

    return Dict([V₀(rg)[i] => round(Bool, vp[i]) + round(Bool, va[i] - vp[i]) * im for i in eachindex(V₀(rg))])
end

function get_polarity_classification_dict(rg::RootGraph, model::JuMP.Model)
    e₊ = value.(model[:e₊])

    return Dict([E(rg)[i] => round(Bool, e₊[i]) for i in eachindex(E(rg))])
end

# multiple graphs

function get_edge_classification_dict(rgs::Vector{RootGraph{T, U}}, models::Vector{<:JuMP.Model}) where {T, U}
    return Dict([i => get_edge_classification_dict(rg, model) for (i, (rg, model)) in enumerate(zip(rgs, models))])
end

function get_segment_classification_dict(rgs::Vector{RootGraph{T, U}}, models::Vector{<:JuMP.Model}) where {T, U}
    return Dict([i => get_segment_classification_dict(rg, model) for (i, (rg, model)) in enumerate(zip(rgs, models))])
end
