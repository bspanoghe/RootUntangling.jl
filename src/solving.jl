function solve_rsa(sg::SuperGraph; optimizer, add_momentum::Bool, time_limit = missing, hotstart_time = 0,
        num_roots::Integer = 1, ρₐ = 0.01, ρₕ = 0.9, ρₘ_max = 0.75, ρₙₙ_max = 0.9, ρᵧ_max = 0.5, ϵ = 1e-2, α_down = -pi/2
    )

    # check for NN prediction data #! remove for final version
    NN_pred = pred_primary(Eₕ₀(sg)[1]) |> !ismissing
    
    # name special vertices
    vₐ = V₊(sg)[1]
    vₑ = V₊(sg)[2] # extinction == disappearance
    vₛ = V₊(sg)[3]

    # define model
    model = Model(optimizer)

    # define the model variables
    n_v = length(V₀(sg))
    n_e = length(E(sg))
    n_he = length(Eₕ₀(sg))
    add_momentum && (n_c = length(E₂(sg)))

    @variable(model, va[1:n_v], Bin) # is vertex active (part of the root)
    @variable(model, vp[1:n_v], Bin) # is vertex part of the primary root

    @variable(model, ea[1:n_e], Bin) # is the edge active (part of the root)
    @variable(model, ep[1:n_e], Bin) # is the edge part of the primary root
    @variable(model, e₊[1:n_e], Bin) # is this a positive edge (should it follow the natural polarity of the edge)
        # note: the natural polarity of an edge is defined as going from the vertex with the lowest id to the one with the highest id

    add_momentum && @variable(model, f[1:n_c], Bin) # are these edges part of the same root

    @variable(model, he[1:n_he], Bin) # is this hyperedge active
    NN_pred && @variable(model, heₚ[1:n_he], Bin) # is this a primary hyperedge

    # connect model variables to graph's edges
    va2f = Dict([V₀(sg)[i] => va[i] for i in eachindex(V₀(sg))])
    vp2f = Dict([V₀(sg)[i] => vp[i] for i in eachindex(V₀(sg))])

    ea2f = Dict([E(sg)[i] => ea[i] for i in eachindex(E(sg))])
    ep2f = Dict([E(sg)[i] => ep[i] for i in eachindex(E(sg))])
    e₊2f = Dict([E(sg)[i] => e₊[i] for i in eachindex(E(sg))])

    add_momentum && (c2f = Dict([E₂(sg)[i] => f[i] for i in eachindex(E₂(sg))]))

    hea2f = Dict([Eₕ₀(sg)[i] => he[i] for i in eachindex(Eₕ₀(sg))])
    NN_pred && (hep2f = Dict(Eₕ₀(sg)[i] => heₚ[i] for i in eachindex(Eₕ₀(sg))))

    # define objective
    @objective(
        model,
        Max,
        # appearance penalties
        sum(
            ea2f[e] * log(ρₐ / (1 - ρₐ))
            for e in E(vₐ)
        ) +
        # standard hyperedges should be active
        sum(
            hea2f[he] * log(ρₕ / (1 - ρₕ))
            for he in Eₕ₀(sg)
        ) + 
        # similar angles
        (
            !add_momentum ? 0 : sum(
            c2f[c] * log(ρₘ(sg, v, c; ρₘ_max, ϵ) / (1 - ρₘ(sg, v, c; ρₘ_max, ϵ)))
            for v in V₀(sg) for c in E₂(v) if !any([is_augmented(e) for e in c])
            )
        ) +
        # gravitropy (needs to be split up into two sums to remain a linear objective)
        sum(
            e₊2f[e] * log(ρᵧ(sg, e, α_down, false; ρᵧ_max, ϵ) / (1 - ρᵧ(sg, e, α_down, false; ρᵧ_max, ϵ)))                
            for e in E₀(sg)
        ) +
        sum(
            (ea2f[e] - e₊2f[e]) * log(ρᵧ(sg, e, α_down, true; ρᵧ_max, ϵ) / (1 - ρᵧ(sg, e, α_down, true; ρᵧ_max, ϵ)))                
            for e in E₀(sg)
        ) +
        # NN pred
        (
            !NN_pred ? 0 : sum(
                hep2f[he] * log(ρₙₙ(he; ρₙₙ_max, ϵ) / (1 - ρₙₙ(he; ρₙₙ_max, ϵ)))
                for he in Eₕ₀(sg)
            )
        ) # +
        # similar predicted thickness
        # (
        #     !add_momentum ? 0 : sum(
        #     c2f[c] * log(ρₜ(sg, c; ρₜ_max, ϵ) / (1 - ρₜ(sg, c; ρₜ_max, ϵ)))
        #     for c in E₂(sg) if !any([is_augmented(e) for e in c])
        #     )
        # )
    )
    
    # # define constraints

    # ## Flow formalism

    # ### For a standard vertex:
    for v in V₀(sg)
        # It can only be classified as primary if active
        @constraint(model, va2f[v] >= vp2f[v])
        # It is an active vertex ⇔ it is connected to two active edges
        @constraint(model, 2*va2f[v] == sum(ea2f[e] for e in E(v)))
        # It is an active primary vertex ⇔ it is connected to two active primary edges
        @constraint(model, 2*vp2f[v] == sum(ep2f[e] for e in E(v)))
        # It has an incoming and an outgoing edge
        @constraint(model, sum((e₊2f[e] - (ea2f[e] - e₊2f[e]))*polarity(e, v) for e in E(v)) == 0)

        # It is active ⇔ It has one active connection
        add_momentum && @constraint(model, va2f[v] == sum(c2f[c] for c in E₂(v)))
    end

    # ### For any edge:
    for e in E(sg)
        # It can only be classified if active
        @constraint(model, ea2f[e] >= ep2f[e])
        @constraint(model, ea2f[e] >= e₊2f[e])
    end

    # ### For a connection:
    add_momentum && for c in E₂(sg)
        # It is active => its edges are active
        @constraint(model, sum(ea2f[e] for e in c) - 1 <= c2f[c])
    end

    # ### For a hyperedge:
    for he in Eₕ₀(sg)
        # It is active => at least one of its edges is active
        @constraint(model, hea2f[he] <= sum(ea2f[e] for e in E(sg, he)))
        # It is classified as primary => at least one of its edges is classified as primary
        NN_pred && @constraint(model, hep2f[he] <= sum(ep2f[e] for e in E(sg, he))) #! model can still choose to not classify he as primary
    end

    # ## Special vertices

    # note: all augmented edges have a positive polarity by design
    # The appearance vertex has no incoming edges (edges are either inactive or follow natural polarity)
    @constraint(model, sum((ea2f[e] - e₊2f[e]) for e in E(vₐ)) == 0)
    # The extinction vertex has no outgoing edges (edges are either inactive or opposite natural polarity)
    @constraint(model, sum(e₊2f[e] for e in E(vₑ)) == 0)
    # The splitting vertex has no incoming edges (edges are either inactive or follow natural polarity)
    @constraint(model, sum((ea2f[e] - e₊2f[e]) for e in E(vₛ)) == 0)
    
    # ## Prerequisite for division

    # A vertex can only split if its hypervertex is part of the primary root 
    # i.e. lateral root segments can only appear in a hypervertex with an edge classified as a primary root
    for v in inner_vertices(sg) # outer nodes can never split
        e_vₛ = edges(v)[findfirst(e -> id(vₛ) ∈ vertices(e), edges(v))] # edge between v and vₛ
        @constraint(model, ea2f[e_vₛ] <= sum(vp2f[v_roommate] for v_roommate in V(sg, Vₕ(v))))
    end

    # Primary root segments cannot form from division
    for e in E(vₛ)
        @constraint(model, ep2f[e] == 0)
    end

    # ## Extras

    if !ismissing(num_roots)
        # The appearance vertex is connected to the primary root with one edge per root
        @constraint(model, sum(ep2f[e] for e in E(vₐ)) == num_roots)
        # The disappearance vertex is connected to the primary root with one edge per root
        @constraint(model, sum(ep2f[e] for e in E(vₑ)) == num_roots)
    end

    # A standard vertex may not have two active edges to the same hypervertex (not required if momentum is added)
    add_momentum || for v in V₀(sg)
        for he in Eₕ₀(v)
            es = filter(e -> id(v) in vertices(e), E(sg, he))
            @constraint(model, sum(ea2f[e] for e in es) <= 1)
        end
    end

    # # Hotstart

    if (hotstart_time > 0) && add_momentum
        @info "Running hotstart"
        start_model = solve_rsa(sg; optimizer, add_momentum = false, time_limit = hotstart_time, num_roots, ρₐ, ρₕ, ϵ)
        vars = all_variables(start_model)
        sols = value.(vars)
        for (start_var, start_sol) in zip(vars, sols)
            varname = JuMP.name(start_var)
            isempty(varname) && continue
            var_current = variable_by_name(model, varname)
            if var_current !== nothing
                set_start_value(var_current, start_sol)
            end
        end
    end
    
    # # extra solver options
    ismissing(time_limit) || set_time_limit_sec(model, time_limit)

    # # solve
    optimize!(model)
    # @assert is_solved_and_feasible(model)

    return model
end

function solve_rsa(sgs::Vector{SuperGraph{T, U}}; kwargs...) where {T, U}
    models = Vector{JuMP.Model}(undef, length(sgs))

    for (i, sg) in enumerate(sgs)
        @info "Solving graph $i/$(length(sgs))"
        models[i] = solve_rsa(sg; kwargs...)
    end
    
    return models
end

# # objective function probabilities
# prevent probabilities from reaching 0 or 1
bound(p; ϵ = 1e-9) = ϵ/2 + (1-ϵ) * p

# change in angle probability
ρₘ(sg, v, c; ρₘ_max, ϵ) = ρₘ_max * angle_dissimilarity(sg, c..., id(v)) |> p -> bound(p; ϵ)

# gravitropic growth probability
ρᵧ(sg, e, α_down, reverse_order; ρᵧ_max, ϵ) = (
    ρᵧ_max * (1 + cosine_similarity(sg, e, α_down; reverse_order)) / 2
) |> p -> bound(p; ϵ)

# NN pred primary probability
ρₙₙ(he; ρₙₙ_max, ϵ) = ρₙₙ_max * pred_primary(he) |> p -> bound(p; ϵ)














# function solve_rsa_quadratic(sg::SuperGraph; optimizer, time_limit = missing, hotstart_time = 0,
#         num_roots::Integer = 1, ρₐ = 0.01, ρₕ = 0.9, ρₘ_max = 0.75, ρₙₙ_max = 0.75, ρᵧ_max = 0.5, ϵ = 1e-2, α_down = -pi/2
#     )

#     # check for NN prediction data #! remove for final version
#     NN_pred = pred_primary(Eₕ₀(sg)[1]) |> !ismissing
    
#     # name special vertices
#     vₐ = V₊(sg)[1]
#     vₑ = V₊(sg)[2] # extinction == disappearance
#     vₛ = V₊(sg)[3]

#     # define model
#     model = Model(optimizer)

#     # define the model variables
#     n_v = length(V₀(sg))
#     n_e = length(E(sg))
#     n_he = length(Eₕ₀(sg))

#     @variable(model, v[1:n_v], Bin) # is vertex active (part of the root)
#     @variable(model, vₚ[1:n_v], Bin) # is vertex part of the primary root
#     @variable(model, vₗ[1:n_v], Bin) # is vertex part of a lateral root

#     @variable(model, e[1:n_e], Bin) # is the edge active (part of the root)
#     @variable(model, eₚ[1:n_e], Bin) # is the edge part of the primary root
#     @variable(model, eₗ[1:n_e], Bin) # is the edge part of a lateral root
#     @variable(model, e₊[1:n_e], Bin) # is this a positive edge (should it follow the natural polarity of the edge)
#     @variable(model, e₋[1:n_e], Bin) # is this a negative edge (should it be opposite the natural polarity of the edge)
#         # note: the natural polarity of an edge is defined as going from the vertex with the lowest id to the one with the highest id

#     @variable(model, he[1:n_he], Bin) # is this hyperedge active
#     NN_pred && @variable(model, heₚ[1:n_he], Bin)

#     # connect model variables to graph's edges
    
#     v2f = Dict([V₀(sg)[i] => v[i] for i in eachindex(V₀(sg))])
#     vp2f = Dict([V₀(sg)[i] => vₚ[i] for i in eachindex(V₀(sg))])
#     vl2f = Dict([V₀(sg)[i] => vₗ[i] for i in eachindex(V₀(sg))])

#     e2f = Dict([E(sg)[i] => e[i] for i in eachindex(E(sg))])
#     ep2f = Dict([E(sg)[i] => eₚ[i] for i in eachindex(E(sg))])
#     el2f = Dict([E(sg)[i] => eₗ[i] for i in eachindex(E(sg))])
#     e₊2f = Dict([E(sg)[i] => e₊[i] for i in eachindex(E(sg))])
#     e₋2f = Dict([E(sg)[i] => e₋[i] for i in eachindex(E(sg))])

#     he2f = Dict([Eₕ₀(sg)[i] => he[i] for i in eachindex(Eₕ₀(sg))])
#     NN_pred && (hep2f = Dict(Eₕ₀(sg)[i] => heₚ[i] for i in eachindex(Eₕ₀(sg))))

#     # define objective
#     @objective(
#         model,
#         Max,
#         # appearance penalties
#         sum(
#             e2f[e] * log(ρₐ / (1 - ρₐ))
#             for e in E(vₐ)
#         ) +
#         # standard hyperedges should be active
#         sum(
#             he2f[he] * log(ρₕ / (1 - ρₕ))
#             for he in Eₕ₀(sg)
#         ) + 
#         # similar angles
#         sum(
#             e2f[e1] * e2f[e2] * log(ρₘ(sg, v, [e1, e2]; ρₘ_max, ϵ) / (1 - ρₘ(sg, v, [e1, e2]; ρₘ_max, ϵ)))
#             for v in V₀(sg) for (e1, e2) in E₂(v) if !any([is_augmented(e) for e in [e1, e2]])
#         ) +
#         # gravitropy (needs to be split up into two sums to remain a linear objective)
#         sum(
#             e₊2f[e] * log(ρᵧ(sg, e, α_down, false; ρᵧ_max, ϵ) / (1 - ρᵧ(sg, e, α_down, false; ρᵧ_max, ϵ)))                
#             for e in E₀(sg)
#         ) +
#         sum(
#             e₋2f[e] * log(ρᵧ(sg, e, α_down, true; ρᵧ_max, ϵ) / (1 - ρᵧ(sg, e, α_down, true; ρᵧ_max, ϵ)))                
#             for e in E₀(sg)
#         ) +
#         # NN pred
#         (
#             !NN_pred ? 0 : sum(
#                 hep2f[he] * log(ρₙₙ(he; ρₙₙ_max, ϵ) / (1 - ρₙₙ(he; ρₙₙ_max, ϵ)))
#                 for he in Eₕ₀(sg)
#             )
#         )
#     )
    
#     # # define constraints

#     # ## Flow formalism

#     # ### For a standard vertex:
#     for v in V₀(sg)
#         # It is active ⇔ it has one active primary/lateral class
#         @constraint(model, v2f[v] == vp2f[v] + vl2f[v])
#         # It is an active primary vertex ⇔ it is connected to two active primary edges
#         @constraint(model, 2*vp2f[v] == sum(ep2f[e] for e in E(v)))
#         # It is an active lateral vertex ⇔ it is connected to two active lateral edges
#         @constraint(model, 2*vl2f[v] == sum(el2f[e] for e in E(v)))
#         # It has an incoming and an outgoing edge
#         @constraint(model, sum((e₊2f[e] - e₋2f[e])*polarity(e, v) for e in E(v)) == 0)
#     end

#     # ### For any edge:
#     for e in E(sg)
#         # It is active ⇔ it has one active primary/lateral class
#         @constraint(model, e2f[e] == ep2f[e] + el2f[e])
#         # It is active ⇔ it has one active outgoing/incoming class
#         @constraint(model, e2f[e] == e₊2f[e] + e₋2f[e])
#     end

#     # ### For a hyperedge:
#     for he in Eₕ₀(sg)
#         # It is active => at least one of its edges is active
#         @constraint(model, he2f[he] <= sum(e2f[e] for e in E(sg, he)))
#         # It is classified as primary => at least one of its edges is classified as primary
#         NN_pred && @constraint(model, hep2f[he] <= sum(ep2f[e] for e in E(sg, he))) #! model can still choose to not classify he as primary
#     end

#     # ## Special vertices

#     # note: all augmented edges have a positive polarity by design
#     # The appearance vertex has no incoming edges (edges are either inactive or follow natural polarity)
#     @constraint(model, sum(e₋2f[e] for e in E(vₐ)) == 0)
#     # The extinction vertex has no outgoing edges (edges are either inactive or opposite natural polarity)
#     @constraint(model, sum(e₊2f[e] for e in E(vₑ)) == 0)
#     # The splitting vertex has no incoming edges (edges are either inactive or follow natural polarity)
#     @constraint(model, sum(e₋2f[e] for e in E(vₛ)) == 0)
    
#     # ## Prerequisite for division

#     # A vertex can only split if its hypervertex is part of the primary root 
#     # i.e. lateral root segments can only appear in a hypervertex with an edge classified as a primary root
#     for v in inner_vertices(sg) # outer nodes can never split
#         e_vₛ = edges(v)[findfirst(e -> id(vₛ) ∈ vertices(e), edges(v))] # edge between v and vₛ
#         @constraint(model, el2f[e_vₛ] <= sum(vp2f[v_roommate] for v_roommate in V(sg, Vₕ(v))))
#     end

#     # Primary root segments cannot form from division
#     for e in E(vₛ)
#         @constraint(model, ep2f[e] == 0)
#     end

#     # ## Exclusion of conflicting hypotheses

#     # For every hypervertex, at most one vertex per exclusion set may be active
#     # i.e. only one hypothesis is true for how many root segments are crossing that hypervertex
#     for v_h in Vₕ₀(sg)
#         for S_l in exclusion_sets(v_h)
#             @constraint(model, sum(v2f[v] for v in V₀(sg)[S_l]) <= 1)
#         end
#     end

#     # ## Extras

#     if !ismissing(num_roots)
#         # The appearance vertex is connected to the primary root with one edge per root
#         @constraint(model, sum(ep2f[e] for e in E(vₐ)) == num_roots)
#         # The disappearance vertex is connected to the primary root with one edge per root
#         @constraint(model, sum(ep2f[e] for e in E(vₑ)) == num_roots)
#     end

#     # # Hotstart

#     if (hotstart_time > 0)
#         @info "Running hotstart"
#         start_model = solve_rsa(sg; optimizer, add_momentum = false, time_limit = hotstart_time, num_roots, ρₐ, ρₕ, ϵ)
#         vars = all_variables(start_model)
#         sols = value.(vars)
#         for (start_var, start_sol) in zip(vars, sols)
#             varname = JuMP.name(start_var)
#             isempty(varname) && continue
#             var_current = variable_by_name(model, varname)
#             if var_current !== nothing
#                 set_start_value(var_current, start_sol)
#             end
#         end
#     end
    
#     # # extra solver options
#     ismissing(time_limit) || set_time_limit_sec(model, time_limit)

#     # # solve
#     optimize!(model)
#     # @assert is_solved_and_feasible(model)

#     return model
# end
