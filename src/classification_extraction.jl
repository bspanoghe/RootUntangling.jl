function get_re_classification_dict(rg::RootGraph, model::JuMP.Model)
    ea = value.(model[:ea])
    ep = value.(model[:ep])

    re_classification_dict = [E(rg)[i] => round(Bool, value(ep[i])) + round(Bool, value(ea[i]) - value(ep[i])) * im for i in eachindex(E(rg))] |> Dict

    return re_classification_dict
end

function get_rv_classification_dict(rg::RootGraph, model::JuMP.Model)
    va = value.(model[:va])
    vp = value.(model[:vp])

    rv_classification_dict = [V₀(rg)[i] => round(Bool, value(vp[i])) + round(Bool, value(va[i]) - value(vp[i])) * im for i in eachindex(V₀(rg))] |> Dict

    return rv_classification_dict
end