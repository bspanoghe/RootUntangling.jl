function get_result_dict(rg::RootGraph, model::JuMP.Model)
    ep = value.(model[:ep])
    el = value.(model[:el])
    ep₊ = value.(model[:ep₊])
    el₊ = value.(model[:el₊])

    return [
        E(rg)[i] => (
            ep = round(Int64, ep[i]),
            el = round(Int64, el[i]),
            ep₊ = round(Int64, ep₊[i]),
            ep₋ = round(Int64, ep[i] - ep₊[i]),
            el₊ = round(Int64, el₊[i]),
            el₋ = round(Int64, el[i] - el₊[i]),
            epa = ep[i] > 0,
            ela = el[i] > 0,
            e₊ = round(Int64, ep₊[i] + el₊[i]),
            e₋ = round(Int64, ep[i] + el[i] - (ep₊[i] + el₊[i])),
            en = round(Int64, ep[i] + el[i]),
        )
        for i in eachindex(E(rg))
    ] |> Dict
end