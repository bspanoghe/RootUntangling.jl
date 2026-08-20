"""
    solve_rsa(
        sg::SuperGraph; optimizer, add_momentum::Bool = true, time_limit = missing, hotstart_time = 0,
        num_roots::Integer = 1, ρₐ = 0.01, ρₕ = 0.95, ρₘ_max = 0.75, ρₙₙ_max = 0.9, ρᵧ_max = 0.5, α_down = -pi/2, ϵ = 1e-5
    )

Solve which root system architecture is represented by the graph `sg`.

The problem is formulated as a Integer Quadratic Program (IQP) written to allow Integer Linear Program (ILP) solvers.

# General keyword arguments
- `optimizer`: The JuMP.jl-compatible ILP optimizer to use.
- `add_momentum`: Minimize angle differences between successive pieces of a root?
- `time_limit`: Time limit of the solver in seconds.
- `hotstart_time`: Time limit of hotstart in seconds. Setting to 0 will disable hotstarting.
- `num_roots`: The amount of root systems present in the graph.
# Solver parameters
- `ρₐ`: The probability of a root appearing without division from the main root.
- `ρₕ`: The probability that an edge truly contains at least one root.
- `ρₘ_max`: The weight given to angle differences, defined as the maximum probability that two succesive edges are the same root if there is no change in angle betweem them.
- `ρₙₙ_max`: The weight given to neural network classfications, defined as the maximum probability that a root is a primary root if classified as such by the neural network.
- `ρᵧ_max`: The weight given to gravitropy, defined as the maximum probability that a downward edge is a root.
- `α_down`: The angle pointing down.
- `ϵ`: The strength of the bound preventing probabilities from reaching 0 or 1 for numerical stability.
"""
function solve_rsa(
        sg::SuperGraph; optimizer, add_momentum::Bool = true, time_limit = missing, hotstart_time = 0,
        num_roots::Integer = 1, ρₐ = 0.01, ρₕ = 0.97, ρₘ_max = 0.75, ρₙₙ_max = 0.9, ρᵧ_max = 0.5, ρₒ_base = exp(1),
        α_down = -pi / 2, ϵ = 1e-3
    )

    @assert ρₒ_base >= 1 "The base for the overlap probability must be greater or equal than 1."

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
        ) +
            # overlap probability
            sum(
            va2f[v] * log(ρ₀(v; ρₒ_base, ϵ) / (1 - ρ₀(v; ρₒ_base, ϵ)))
                for v in V₀(sg)
        )
    )

    # # define constraints

    # ## Flow formalism

    # ### For a standard vertex:
    for v in V₀(sg)
        # It can only be classified as primary if active
        @constraint(model, va2f[v] >= vp2f[v])
        # It is an active vertex ⇔ it is connected to two active edges
        @constraint(model, 2 * va2f[v] == sum(ea2f[e] for e in E(v)))
        # It is an active primary vertex ⇔ it is connected to two active primary edges
        @constraint(model, 2 * vp2f[v] == sum(ep2f[e] for e in E(v)))
        # It has an incoming and an outgoing edge
        @constraint(model, sum((e₊2f[e] - (ea2f[e] - e₊2f[e])) * polarity(e, v) for e in E(v)) == 0)

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
    add_momentum || for v in V₀(sg), he in Eₕ₀(v)
        es = filter(e -> id(v) in vertices(e), E(sg, he))
        @constraint(model, sum(ea2f[e] for e in es) <= 1)
    end

    # # Hotstart

    if (hotstart_time > 0) && add_momentum
        @info "Running hotstart"
        start_model = solve_rsa(sg; optimizer, add_momentum = false, time_limit = hotstart_time,
            num_roots, ρₐ, ρₕ, ρₘ_max, ρₙₙ_max, ρᵧ_max, ρₒ_base, α_down, ϵ)

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
bound(p; ϵ = 1.0e-9) = ϵ / 2 + (1 - ϵ) * p

# change in angle probability
ρₘ(sg, v, c; ρₘ_max, ϵ) = ρₘ_max * angle_dissimilarity(sg, c..., id(v)) |> p -> bound(p; ϵ)

# gravitropic growth probability
ρᵧ(sg, e, α_down, reverse_order; ρᵧ_max, ϵ) = (
    ρᵧ_max * (1 + cosine_similarity(sg, e, α_down; reverse_order)) / 2
) |> p -> bound(p; ϵ)

# NN pred primary probability
ρₙₙ(he; ρₙₙ_max, ϵ) = ρₙₙ_max * pred_primary(he) |> p -> bound(p; ϵ)

# root overlap probability
ρ₀(sv::SingularVertex; ρₒ_base, ϵ) = 0.5 * ρₒ_base^(-order(sv)) |> p -> bound(p; ϵ)
order(sv::SingularVertex) = id(sv) - vertices(hypervertex(sv))[1]
