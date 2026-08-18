# # classification of graph edges/vertices
# one graph

function get_se_classification_dict(sg::SuperGraph, model::JuMP.Model)
    ea = [var for var in all_variables(model) if !isnothing(match(r"^ea\[\d+\]$", JuMP.name(var)))]
    ep = [var for var in all_variables(model) if !isnothing(match(r"^ep\[\d+\]$", JuMP.name(var)))]

    se_classification_dict = [E(sg)[i] => round(Bool, value(ep[i])) + round(Bool, value(ea[i]) - value(ep[i]))*im for i in eachindex(E(sg))] |> Dict

    return se_classification_dict
end

function get_he_classification_dict(sg::SuperGraph, model::JuMP.Model)
    se_classification_dict = get_se_classification_dict(sg, model)
    he_classification_dict = [he => sum([se_classification_dict[e] for e in E(sg, he)]) for he in Eₕ(sg)] |> Dict

    return he_classification_dict
end

function get_sv_classification_dict(sg::SuperGraph, model::JuMP.Model)
    va = [var for var in all_variables(model) if !isnothing(match(r"^va\[\d+\]$", JuMP.name(var)))]
    vp = [var for var in all_variables(model) if !isnothing(match(r"^vp\[\d+\]$", JuMP.name(var)))]

    sv_classification_dict = [V₀(sg)[i] => round(Bool, value(vp[i])) + round(Bool, value(va[i]) - value(vp[i]))*im for i in eachindex(V₀(sg))] |> Dict

    return sv_classification_dict
end

function get_hv_classification_dict(sg::SuperGraph, model::JuMP.Model)
    sv_classification_dict = get_sv_classification_dict(sg, model)
    hv_classification_dict = [hv => sum([sv_classification_dict[v] for v in V(sg, hv)]) for hv in Vₕ₀(sg)] |> Dict

    return hv_classification_dict
end

function get_polarity_classification_dict(sg::SuperGraph, model::JuMP.Model)
    e₊ = [var for var in all_variables(model) if !isnothing(match(r"^e₊\[\d+\]$", JuMP.name(var)))]
    
    polarity_classification_dict = [E(sg)[i] => round(Bool, value(e₊[i])) for i in eachindex(E(sg))] |> Dict

    return polarity_classification_dict
end

# multiple graphs

function get_se_classification_dict(sgs::Vector{SuperGraph{T, U}}, models::Vector{<:JuMP.Model}; kwargs...) where {T, U}
    classification_dict = Dict{Int64, Dict{SingularEdge{T, U}, Complex{Int64}}}()

    for (i, (sg, model)) in enumerate(zip(sgs, models))
        classification_dict[i] = get_se_classification_dict(sg, model)
    end
    
    return classification_dict
end

function get_he_classification_dict(sgs::Vector{SuperGraph{T, U}}, models::Vector{<:JuMP.Model}; kwargs...) where {T, U}
    classification_dict = Dict{Int64, Dict{HyperEdge{T, U}, Complex{Int64}}}()

    for (i, (sg, model)) in enumerate(zip(sgs, models))
        classification_dict[i] = get_he_classification_dict(sg, model)
    end
    
    return classification_dict
end