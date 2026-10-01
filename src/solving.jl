"""
    solve_rsa(
        rg::RootGraph; optimizer, time_limit = missing, hotstart_time = 0,
        num_roots::Integer = 1, ρₐ = 0.01, ρₕ = 0.95, ρₘ_max = 0.75, ρₙₙ_max = 0.9, ρᵧ_max = 0.5, α_down = -pi/2, ϵ = 1e-5
    )

Solve which root system architecture is represented by the graph `rg`.

The problem is formulated as a Integer Quadratic Program (IQP) written to allow Integer Linear Program (ILP) solvers.

# General keyword arguments
- `optimizer`: The JuMP.jl-compatible ILP optimizer to use.
- `time_limit`: Time limit of the solver in seconds.
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
        rg::RootGraph; optimizer, time_limit = missing, max_overlapping = 3,
        num_roots::Integer = 1, ρₐ = 0.01, ρₕ = 0.97, ρₒ = 0.1, ρₘ_max = 0.75, ρₙₙ_max = 0.9, ρᵧ_max = 0.5,
        w_gp = 1.0, w_gl = 1.0, α_down = -pi / 2, ϵ = 1e-3
    )

    # check for NN prediction data #! remove for final version
    NN_pred = pred_primary(E₀(rg)[1]) |> !ismissing

    # name special vertices
    vₐ = V₊(rg)[1]
    vₑ = V₊(rg)[2] # extinction == disappearance
    vₛ = V₊(rg)[3]

    # define model
    model = Model(optimizer)

    # define the model variables
    connections = E₂(rg)

    n_e = length(E(rg))
    n_c = length(connections)

    @variable(model, ea[1:n_e], Bin) # does edge contain roots
    @variable(model, epa[1:n_e], Bin) # does edge contain primary roots
    @variable(model, ep₊[1:n_e], Int, lower_bound = 0, upper_bound = max_overlapping) # number of roots in edge following positive direction
    @variable(model, ep₋[1:n_e], Int, lower_bound = 0, upper_bound = max_overlapping) # number of primary roots in edge
    @variable(model, el₊[1:n_e], Int, lower_bound = 0, upper_bound = max_overlapping) # number of roots in edge following positive direction
    @variable(model, el₋[1:n_e], Int, lower_bound = 0, upper_bound = max_overlapping) # number of lateral roots in edge

    @variable(model, cp₊[1:n_c], Int, lower_bound = 0, upper_bound = max_overlapping) # number of edges in primary connection
    @variable(model, cp₋[1:n_c], Int, lower_bound = 0, upper_bound = max_overlapping) # number of edges in primary connection
    @variable(model, cl₊[1:n_c], Int, lower_bound = 0, upper_bound = max_overlapping) # number of edges in lateral connection
    @variable(model, cl₋[1:n_c], Int, lower_bound = 0, upper_bound = max_overlapping) # number of edges in lateral connection

    # connect model variables to graph's edges
    ea2f = Dict([E(rg)[i] => ea[i] for i in eachindex(E(rg))])
    epa2f = Dict([E(rg)[i] => epa[i] for i in eachindex(E(rg))])
    ep₊2f = Dict([E(rg)[i] => ep₊[i] for i in eachindex(E(rg))])
    ep₋2f = Dict([E(rg)[i] => ep₋[i] for i in eachindex(E(rg))])
    el₊2f = Dict([E(rg)[i] => el₊[i] for i in eachindex(E(rg))])
    el₋2f = Dict([E(rg)[i] => el₋[i] for i in eachindex(E(rg))])

    cp₊2f = Dict([connections[i] => cp₊[i] for i in eachindex(connections)])
    cp₋2f = Dict([connections[i] => cp₋[i] for i in eachindex(connections)])
    cl₊2f = Dict([connections[i] => cl₊[i] for i in eachindex(connections)])
    cl₋2f = Dict([connections[i] => cl₋[i] for i in eachindex(connections)])

    # define objective
    @objective(
        model,
        Max,
        # appearance penalties
        sum(
            (ep₊2f[e] + ep₋2f[e] + el₊2f[e] + el₋2f[e]) * log(ρₐ / (1 - ρₐ))
            for e in E(vₐ)
        ) +
        # standard rootedges should be active
        sum(
            ea2f[e] * log(ρₕ / (1 - ρₕ))
            for e in E₀(rg)
        ) + 
        # overlap probability
        sum(
            (ep₊2f[e] + ep₋2f[e] + el₊2f[e] + el₋2f[e]) * log(ρₒ / (1 - ρₒ)) #! -1?
            for e in E₀(rg)
        ) +
        # primary gravitropy (needs to be split up into two sums to remain a linear objective)
        w_gp * sum(
            ep₊2f[e] * log(ρᵧ(rg, e, α_down, false; ρᵧ_max, ϵ) / (1 - ρᵧ(rg, e, α_down, false; ρᵧ_max, ϵ)))
            for e in E₀(rg)
        ) +
        w_gp * sum(
            ep₋2f[e] * log(ρᵧ(rg, e, α_down, true; ρᵧ_max, ϵ) / (1 - ρᵧ(rg, e, α_down, true; ρᵧ_max, ϵ)))
            for e in E₀(rg)
        ) + 
        # lateral gravitropy (needs to be split up into two sums to remain a linear objective)
        w_gl * sum(
            el₊2f[e] * log(ρᵧ(rg, e, α_down, false; ρᵧ_max, ϵ) / (1 - ρᵧ(rg, e, α_down, false; ρᵧ_max, ϵ)))
            for e in E₀(rg)
        ) +
        w_gl * sum(
            el₋2f[e] * log(ρᵧ(rg, e, α_down, true; ρᵧ_max, ϵ) / (1 - ρᵧ(rg, e, α_down, true; ρᵧ_max, ϵ)))
            for e in E₀(rg)
        ) + 
        # angle differences
        sum(
            (cp₊2f[c] + cp₋2f[c] + cl₊2f[c] + cl₋2f[c]) * log(ρₘ(rg, v, c; ρₘ_max, ϵ) / (1 - ρₘ(rg, v, c; ρₘ_max, ϵ)))
            for v in V₀(rg) for c in E₂(v) if !any([is_augmented(e) for e in c])
        )
    )

    # NN pred
    if NN_pred
        set_objective_function(
            model,
            objective_function(model) + sum(
                epa2f[e] * log(ρₙₙ(e; ρₙₙ_max, ϵ) / (1 - ρₙₙ(e; ρₙₙ_max, ϵ)))
                for e in E₀(rg)
            )
        )
    end
    
    # # define constraints

    # ## Flow formalism

    # ### For a standard vertex:
    for v in V₀(rg)
        # For all of its edges:
        for e in E(v)
            # it contains as many roots as the sum of its connections
            @constraint(model, ep₊2f[e] + ep₋2f[e] == sum(cp₊2f[c] + cp₋2f[c] for c in E₂(v, e)))
            @constraint(model, el₊2f[e] + el₋2f[e] == sum(cl₊2f[c] + cl₋2f[c] for c in E₂(v, e)))

            # it has as many roots in a given direction as the sum of its connected roots in that direction
            @constraint(
                model, 
                (ep₊2f[e] - ep₋2f[e]) * direction(v, e) == 
                    sum((cp₊2f[c] - cp₋2f[c]) * direction(c, e) for c in E₂(v, e))
            )
            @constraint(
                model,
                (el₊2f[e] - el₋2f[e]) * direction(v, e) == 
                    sum((cl₊2f[c] - cl₋2f[c]) * direction(c, e) for c in E₂(v, e))
            )
        end
    end

    # ### For any edge:
    for e in E(rg)
        # It can only be classified as active if it contains roots
        @constraint(model, ea2f[e] <= (ep₊2f[e] + ep₋2f[e] + el₊2f[e] + el₋2f[e]))
        @constraint(model, epa2f[e] <= ep₊2f[e] + ep₋2f[e])
    end

    # ## Special vertices

    # The appearance vertex has no incoming edges (edges are either inactive or follow natural polarity)
    @constraint(model, sum(ep₋2f[e] for e in E(vₐ)) == 0)
    @constraint(model, sum(el₋2f[e] for e in E(vₐ)) == 0)
    # The extinction vertex has no outgoing edges (edges are either inactive or opposite natural polarity)
    @constraint(model, sum(ep₊2f[e] for e in E(vₑ)) == 0)
    @constraint(model, sum(el₊2f[e] for e in E(vₑ)) == 0)
    # The splitting vertex has no incoming edges (edges are either inactive or follow natural polarity)
    @constraint(model, sum(ep₋2f[e] for e in E(vₛ)) == 0)
    @constraint(model, sum(el₋2f[e] for e in E(vₛ)) == 0)

    # ## Prerequisite for division

    # A vertex can only split if it's part of the primary root
    for v in inner_vertices(rg) # outer nodes can never split
        e_vₛ = edges(v)[findfirst(e -> id(vₛ) ∈ vertices(e), edges(v))] # edge between v and vₛ
        @constraint(model, el₊2f[e_vₛ] <= sum(ep₊2f[e] + ep₋2f[e] for e in E(v))) #! only allows one split event in vertex (otherwise multiply right side by constant)
    end

    # Primary root segments cannot form from division
    for e in E(vₛ)
        @constraint(model, ep₊2f[e] == 0)
    end

    # ## Extras

    if !ismissing(num_roots)
        # The appearance vertex is connected to the primary root with one edge per root
        @constraint(model, sum(ep₊2f[e] for e in E(vₐ)) == num_roots)
        # The disappearance vertex is connected to the primary root with one edge per root
        @constraint(model, sum(ep₋2f[e] for e in E(vₑ)) == num_roots)
    end

    # # extra solver options
    ismissing(time_limit) || set_time_limit_sec(model, time_limit)

    # # solve
    optimize!(model)
    # @assert is_solved_and_feasible(model)

    return model
end

function solve_rsa(rgs::Vector{RootGraph{T, U}}; kwargs...) where {T, U}
    models = Vector{JuMP.Model}(undef, length(rgs))

    for (i, rg) in enumerate(rgs)
        @info "Solving graph $i/$(length(rgs))"
        models[i] = solve_rsa(rg; kwargs...)
    end

    return models
end

# # objective function probabilities
# prevent probabilities from reaching 0 or 1
bound(p; ϵ = 1.0e-9) = ϵ / 2 + (1 - ϵ) * p

# change in angle probability
ρₘ(rg, v, c; ρₘ_max, ϵ) = ρₘ_max * angle_dissimilarity(rg, c..., id(v)) |> p -> bound(p; ϵ)

# gravitropic growth probability
ρᵧ(rg, e, α_down, reverse_order; ρᵧ_max, ϵ) = (
    ρᵧ_max * (1 + cosine_similarity(rg, e, α_down; reverse_order)) / 2
) |> p -> bound(p; ϵ)

# NN pred primary probability
ρₙₙ(re; ρₙₙ_max, ϵ) = ρₙₙ_max * pred_primary(re) |> p -> bound(p; ϵ)