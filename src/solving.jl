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
        rg::RootGraph; optimizer, time_limit = missing,
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

    n_e = length(E(rg))
    n_c = length(connections)

    #! upper bounds
    @variable(model, ea[1:n_e], Bin) # is edge active
    @variable(model, ep[1:n_e], Int, lower_bound = 0, upper_bound = 3) # number of primary roots in edge
    @variable(model, ep₊[1:n_e], Int, lower_bound = 0, upper_bound = 3) # number of roots in edge following positive direction
    @variable(model, el[1:n_e], Int, lower_bound = 0, upper_bound = 3) # number of lateral roots in edge
    @variable(model, el₊[1:n_e], Int, lower_bound = 0, upper_bound = 3) # number of roots in edge following positive direction

    @variable(model, fp[1:n_c], Int, lower_bound = 0, upper_bound = 3) # number of edges in primary connection
    @variable(model, fl[1:n_c], Int, lower_bound = 0, upper_bound = 3) # number of edges in lateral connection

    # connect model variables to graph's edges
    ea2f = Dict([E(rg)[i] => ea[i] for i in eachindex(E(rg))])
    ep2f = Dict([E(rg)[i] => ep[i] for i in eachindex(E(rg))])
    ep₊2f = Dict([E(rg)[i] => ep₊[i] for i in eachindex(E(rg))])
    el2f = Dict([E(rg)[i] => el[i] for i in eachindex(E(rg))])
    el₊2f = Dict([E(rg)[i] => el₊[i] for i in eachindex(E(rg))])

    cp2f = Dict([connections[i] => fp[i] for i in eachindex(connections)])
    cl2f = Dict([connections[i] => fl[i] for i in eachindex(connections)])

    # define objective
    @objective(
        model,
        Max,
        # appearance penalties
        sum(
            (ep2f[e] + el2f[e]) * log(ρₐ / (1 - ρₐ))
            for e in E(vₐ)
        ) +
        # standard rootedges should be active
        sum(
            ea2f[e] * log(ρₕ / (1 - ρₕ))
            for e in E₀(rg)
        ) + 
        # but not TOO many #!
        -0.4 * sum(
            (ep2f[e] + el2f[e]) * log(ρₕ / (1 - ρₕ)) #!
            for e in E₀(rg)
        ) +
        # gravitropy (needs to be split up into two sums to remain a linear objective)
        sum(
            (ep₊2f[e] + el₊2f[e]) * log(ρᵧ(rg, e, α_down, false; ρᵧ_max, ϵ) / (1 - ρᵧ(rg, e, α_down, false; ρᵧ_max, ϵ)))
            for e in E₀(rg)
        ) +
        sum(
            ((ep2f[e] + el2f[e]) - (ep₊2f[e] + el₊2f[e])) * log(ρᵧ(rg, e, α_down, true; ρᵧ_max, ϵ) / (1 - ρᵧ(rg, e, α_down, true; ρᵧ_max, ϵ)))
            for e in E₀(rg)
        ) + 
        # angle differences
        sum(
            (cp2f[c] + cl2f[c]) * log(ρₘ(rg, v, c; ρₘ_max, ϵ) / (1 - ρₘ(rg, v, c; ρₘ_max, ϵ)))
            for v in V₀(rg) for c in E₂(v) if !any([is_augmented(e) for e in c])
        )
    )

    # NN pred
    if NN_pred
        set_objective_function(
            model,
            objective_function(model) + sum(
                ep2f[e] * log(ρₙₙ(e; ρₙₙ_max, ϵ) / (1 - ρₙₙ(e; ρₙₙ_max, ϵ))) #! stacks w num of primary roots
                for e in E₀(rg)
            )
        )
    end
    
    # # define constraints

    # ## Flow formalism

    # ### For a standard vertex:
    for v in V₀(rg)
        # For all of its standard edges:
        for e in E(v)
            # it contains as many roots as the sum of its connections
            @constraint(model, ep2f[e] == sum(cp2f[c] for c in E₂(v, e)))
            @constraint(model, el2f[e] == sum(cl2f[c] for c in E₂(v, e)))
        end
        # It has equal incoming and outgoing roots
        @constraint(model, sum((ep₊2f[e] - (ep2f[e] - ep₊2f[e])) * direction(v, e) for e in E(v)) == 0)
        @constraint(model, sum((el₊2f[e] - (el2f[e] - el₊2f[e])) * direction(v, e) for e in E(v)) == 0)
    end

    # ### For any edge:
    for e in E(rg)
        # It can only be classified as active if it contains roots
        @constraint(model, ea2f[e] <= (ep2f[e] + el2f[e]))

        # The amount of roots in positive direction is no larger than the amount of total roots
        @constraint(model, ep₊2f[e] <= ep2f[e])
        @constraint(model, el₊2f[e] <= el2f[e])
    end

    # ## Special vertices

    # The appearance vertex has no incoming edges (edges are either inactive or follow natural polarity)
    @constraint(model, sum((ep2f[e] - ep₊2f[e]) for e in E(vₐ)) == 0)
    @constraint(model, sum((el2f[e] - el₊2f[e]) for e in E(vₐ)) == 0)
    # The extinction vertex has no outgoing edges (edges are either inactive or opposite natural polarity)
    @constraint(model, sum(ep₊2f[e] for e in E(vₑ)) == 0)
    @constraint(model, sum(el₊2f[e] for e in E(vₑ)) == 0)
    # The splitting vertex has no incoming edges (edges are either inactive or follow natural polarity)
    @constraint(model, sum((ep2f[e] - ep₊2f[e]) for e in E(vₛ)) == 0)
    @constraint(model, sum((el2f[e] - el₊2f[e]) for e in E(vₛ)) == 0)

    # ## Prerequisite for division

    # A vertex can only split if it's part of the primary root
    for v in inner_vertices(rg) # outer nodes can never split
        e_vₛ = edges(v)[findfirst(e -> id(vₛ) ∈ vertices(e), edges(v))] # edge between v and vₛ
        @constraint(model, el2f[e_vₛ] <= sum(ep2f[e] for e in E(v))) #! only allows one split event in vertex (otherwise multiply right side by constant)
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
