"""
    solve_rsa(
        rg::RootGraph; optimizer, add_momentum::Bool = true, time_limit = missing, hotstart_time = 0,
        num_roots::Integer = 1, ρₐ = 0.01, ρₕ = 0.95, ρₘ_max = 0.75, ρₙₙ_max = 0.9, ρᵧ_max = 0.5, α_down = -pi/2, ϵ = 1e-5
    )

Solve which root system architecture is represented by the graph `rg`.

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
        rg::RootGraph; optimizer, add_momentum::Bool = true, time_limit = missing, hotstart_time = 0,
        num_roots::Integer = 1, ρₐ = 0.01, ρₕ = 0.97, ρₘ_max = 0.75, ρₙₙ_max = 0.9, ρᵧ_max = 0.5, ρₒ_base = exp(1),
        α_down = -pi / 2, ϵ = 1e-3
    )

    @assert ρₒ_base >= 1 "The base for the overlap probability must be greater or equal than 1."

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

    n_v = length(V₀(rg))
    n_e = length(E(rg))
    add_momentum && (n_c = length(connections))

    @variable(model, vp[1:n_v], Bin) # is vertex part of the primary root
    @variable(model, vn[1:n_v], Int, lower_bound = 0, upper_bound = 3) # number of roots in vertex #! ub

    @variable(model, ea[1:n_e], Bin) # is edge active
    @variable(model, ep[1:n_e], Bin) # is the edge part of the primary root
    @variable(model, en[1:n_e], Int, lower_bound = 0, upper_bound = 3) # number of roots in edge #! ub
    @variable(model, e₊[1:n_e], Int, lower_bound = 0, upper_bound = 3) # number of roots in edge following positive direction #! ub

    add_momentum && @variable(model, f[1:n_c], Int, lower_bound = 0, upper_bound = 3) # number of edges in connection #! ub

    # connect model variables to graph's edges
    vn2f = Dict([V₀(rg)[i] => vn[i] for i in eachindex(V₀(rg))])
    vp2f = Dict([V₀(rg)[i] => vp[i] for i in eachindex(V₀(rg))])

    ea2f = Dict([E(rg)[i] => ea[i] for i in eachindex(E(rg))])
    ep2f = Dict([E(rg)[i] => ep[i] for i in eachindex(E(rg))])
    en2f = Dict([E(rg)[i] => en[i] for i in eachindex(E(rg))])
    e₊2f = Dict([E(rg)[i] => e₊[i] for i in eachindex(E(rg))])

    add_momentum && (c2f = Dict([connections[i] => f[i] for i in eachindex(connections)]))

    # define objective
    @objective(
        model,
        Max,
        # appearance penalties
        sum(
            en2f[e] * log(ρₐ / (1 - ρₐ))
            for e in E(vₐ)
        ) +
        # standard rootedges should be active
        sum(
            ea2f[e] * log(ρₕ / (1 - ρₕ))
            for e in E₀(rg)
        ) + 
        # but not TOO many #!
        -0.75 * sum(
            en2f[e] * log(ρₕ / (1 - ρₕ)) #!
            for e in E₀(rg)
        ) +
        # gravitropy (needs to be split up into two sums to remain a linear objective)
        sum(
            e₊2f[e] * log(ρᵧ(rg, e, α_down, false; ρᵧ_max, ϵ) / (1 - ρᵧ(rg, e, α_down, false; ρᵧ_max, ϵ)))
            for e in E₀(rg)
        ) +
        sum(
            (en2f[e] - e₊2f[e]) * log(ρᵧ(rg, e, α_down, true; ρᵧ_max, ϵ) / (1 - ρᵧ(rg, e, α_down, true; ρᵧ_max, ϵ)))
            for e in E₀(rg)
        )
    )

    # similar angles
    if add_momentum
        set_objective_function(
            model,
            objective_function(model) + sum(
                c2f[c] * log(ρₘ(rg, v, c; ρₘ_max, ϵ) / (1 - ρₘ(rg, v, c; ρₘ_max, ϵ)))
                for v in V₀(rg) for c in E₂(v) if !any([is_augmented(e) for e in c])
            )
        )
    end

    # NN pred
    if NN_pred
        set_objective_function(
            model,
            objective_function(model) + sum(
                ep2f[e] * log(ρₙₙ(e; ρₙₙ_max, ϵ) / (1 - ρₙₙ(e; ρₙₙ_max, ϵ)))
                for e in E₀(rg)
            )
        )
    end
    
    # # define constraints

    # ## Flow formalism

    # ### For a standard vertex:
    for v in V₀(rg)
        # It can only be classified as primary if it contains roots
        @constraint(model, vp2f[v] <= vn2f[v])
        # it contains n roots ⇔ it is connected to edges summing to 2n roots (n incoming + n outgoing)
        @constraint(model, 2 * vn2f[v] == sum(en2f[e] for e in E(v)))
        # It is an active primary vertex ⇔ it is connected to two active primary edges
        @constraint(model, 2 * vp2f[v] == sum(ep2f[e] for e in E(v)))
        # It has equal incoming and outgoing roots
        @constraint(model, sum((e₊2f[e] - (en2f[e] - e₊2f[e])) * direction(v, e) for e in E(v)) == 0)
        # The amount of roots passing through equals those in the sum of its connections
        add_momentum && @constraint(model, vn2f[v] == sum(c2f[c] for c in E₂(v)))

        # For all of its standard edges:
        for e in E₀(v)
            # it contains as many roots as the sum of its connections
            add_momentum && @constraint(model, en2f[e] == sum(c2f[c] for c in E₂(v, e)))
        end
    end

    # ### For any edge:
    for e in E(rg)
        # It can only be classified as active if it contains roots
        @constraint(model, ea2f[e] <= en2f[e])
        # It can only be classified as primary if it contains roots
        @constraint(model, ep2f[e] <= en2f[e])
        # The amount of roots in positive direction is no larger than the amount of total roots
        @constraint(model, e₊2f[e] <= en2f[e])
    end

    # ### For a connection:
    add_momentum && for c in connections
        # It is active => its edges are active
        # @constraint(model, sum(ea2f[e] for e in c) - 1 <= c2f[c])
        # It is active => Its edges have the same classification
        # @constraint(model, ep2f[c[1]] == ep2f[c[2]])
    end

    # ## Special vertices

    # The appearance vertex has no incoming edges (edges are either inactive or follow natural polarity)
    @constraint(model, sum((en2f[e] - e₊2f[e]) for e in E(vₐ)) == 0)
    # The extinction vertex has no outgoing edges (edges are either inactive or opposite natural polarity)
    @constraint(model, sum(e₊2f[e] for e in E(vₑ)) == 0)
    # The splitting vertex has no incoming edges (edges are either inactive or follow natural polarity)
    @constraint(model, sum((en2f[e] - e₊2f[e]) for e in E(vₛ)) == 0)

    # ## Prerequisite for division

    # A vertex can only split if it's part of the primary root
    for v in inner_vertices(rg) # outer nodes can never split
        e_vₛ = edges(v)[findfirst(e -> id(vₛ) ∈ vertices(e), edges(v))] # edge between v and vₛ
        @constraint(model, en2f[e_vₛ] <= vp2f[v]) #! only allows one split event in vertex (otherwise multiply right side by constant)
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

    # # Hotstart

    if (hotstart_time > 0) && add_momentum
        @info "Running hotstart"
        start_model = solve_rsa(rg; optimizer, add_momentum = false, time_limit = hotstart_time,
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

# root overlap probability
ρ₀(rv::RootVertex; ρₒ_base, ϵ) = 0.5 * ρₒ_base^(-order(rv)) |> p -> bound(p; ϵ)
order(rv::RootVertex) = id(rv) - vertices(rootvertex(rv))[1]
