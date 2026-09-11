function get_re_classification_dict(rg::RootGraph, model::JuMP.Model)
    en = value.(model[:en])
    ep = value.(model[:ep])

    re_classification_dict = [E(rg)[i] => round(Bool, ep[i]) + (en[i] > 0 && ep[i] < 1) * im for i in eachindex(E(rg))] |> Dict

    return re_classification_dict
end

function get_en_dict(rg::RootGraph, model::JuMP.Model)
    en = value.(model[:en])
    return [E(rg)[i] => round(Int64, en[i]) for i in eachindex(E(rg))] |> Dict
end

function get_polarity_dict(rg::RootGraph, model::JuMP.Model)
    e₊ = value.(model[:e₊])
    return [E(rg)[i] => round(Int64, e₊[i]) for i in eachindex(E(rg))] |> Dict
end