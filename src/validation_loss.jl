function get_loss(rs::RootSystem, vrs::ValRootSystem, root_mapping::Vector{<:Integer}, metric::Function;
        verbose::Bool = false
    )
    #! prevent deepcopy here for efficiency
    rs, vrs = deepcopy(rs), deepcopy(vrs)
    add_dummy_roots!(rs, vrs)

    # primary root
    primary_dist = get_loss(primary(rs), primary(vrs), metric)

    # laterals
    lateral_dist = 0
    for i in eachindex(root_mapping)
        lateral_dist += get_loss(laterals(rs), laterals(vrs), i, root_mapping[i], metric)
        verbose && @info "Loss in root $(i): $(get_loss(laterals(rs), laterals(vrs), i, root_mapping[i], metric))"
    end

    return primary_dist + lateral_dist
end

get_loss(rs::Vector{<:Root}, rvs::Vector{<:ValRoot}, i::Integer, j::Integer, metric::Function) = get_loss(rs[i], rvs[j], metric)

function get_loss(r::Root, rv::ValRoot, metric::Function)
    return length(full_setdiff(metric(r), metric(rv)))
end

full_setdiff(x1::AbstractArray, x2::AbstractArray) = [setdiff(x1, x2); setdiff(x2, x1)]

function get_optimal_mapping(rs::RootSystem, vrs::ValRootSystem, metric::Function; optimizer)
    rs, vrs = deepcopy(rs), deepcopy(vrs)
    add_dummy_roots!(rs, vrs)
    return get_optimal_mapping(laterals(rs), laterals(vrs), metric; optimizer)
end

function get_optimal_mapping(rs::Vector{<:Root}, vrs::Vector{<:ValRoot}, metric::Function; optimizer)
    n = length(rs)

    loss_mat = [
        get_loss(rs, vrs, i, j, metric)
        for i in 1:n, j in 1:n
    ]
  
    return get_optimal_mapping(loss_mat; optimizer)
end

function get_optimal_mapping(loss_mat::Matrix{<:Number}; optimizer)
    n = size(loss_mat, 1)
    model = Model(optimizer)
    @variable(model, mm[1:n, 1:n], Bin)
    @objective(model, Min, sum(loss_mat .* mm))
    for i in 1:n
        @constraint(model, sum(mm[i,j] for j in 1:n) == 1)
    end
    for j in 1:n
        @constraint(model, sum(mm[i,j] for i in 1:n) == 1)
    end

    set_silent(model)
    optimize!(model)
    is_solved_and_feasible(model)
    mm = round.(Bool, value.(model[:mm]))
    return [findfirst(c) for c in eachrow(mm)]
end

function get_optimal_mappings(rss::Vector{<:RootSystem}, vrss::Vector{<:ValRootSystem}, metric::Function; optimizer)
    n_og = length(rss)
    rss, vrss = deepcopy(rss), deepcopy(vrss)
    add_dummy_rootsystems!(rss, vrss)
    n = length(rss)

    mapping_mat = [
        get_optimal_mapping(rss[i], vrss[j], metric; optimizer)
        for i in 1:n, j in 1:n
    ]
    loss_mat = [
        get_loss(rss[i], vrss[j], mapping_mat[i, j], metric)
        for i in 1:n, j in 1:n
    ]

    rs_mapping = get_optimal_mapping(loss_mat; optimizer)
    return rs_mapping[1:n_og], [mapping_mat[i, j] for (i, j) in zip(1:n_og, rs_mapping[1:n_og])]
end

add_dummy_roots!(rs::RootSystem, vrs::ValRootSystem) = add_dummy_roots!(laterals(rs), laterals(vrs))

function add_dummy_roots!(rs::Vector{<:Root}, vrs::Vector{<:ValRoot})
    n_diff = length(rs) - length(vrs)
    n_diff < 0 && add_dummy_roots!(rs, abs(n_diff))
    n_diff > 0 && add_dummy_roots!(vrs, n_diff)

    return nothing
end

function add_dummy_roots!(roots::Vector{Union{CompositeRoot{T, U}, SimpleRoot{T, U}}}, n_diff::Integer) where {T, U}
    append!(roots, fill(SimpleRoot(false, RootFragment(false, RootArc{T, U}[])), n_diff))

    return nothing
end

function add_dummy_roots!(roots::Vector{ValRoot{T, U}}, n_diff::Integer) where {T, U}
    append!(roots, fill(ValRoot(false, RootEdge{T, U}[], U[], U[]), n_diff))

    return nothing
end

function add_dummy_rootsystems!(rss::Vector{<:RootSystem}, vrss::Vector{<:ValRootSystem})
    n_diff = length(rss) - length(vrss)
    n_diff < 0 && add_dummy_rootsystems!(rss, abs(n_diff))
    n_diff > 0 && add_dummy_rootsystems!(vrss, n_diff)

    return nothing
end

function add_dummy_rootsystems!(rss::Vector{RootSystem{T, U}}, n_diff::Integer) where {T, U}
    append!(
        rss, 
        fill(
            RootSystem(
                SimpleRoot(true, RootFragment(false, RootArc{T, U}[])),
                Union{SimpleRoot{T, U}, CompositeRoot{T, U}}[], 
            ),
            n_diff
        )
    )

    return nothing
end

function add_dummy_rootsystems!(vrss::Vector{ValRootSystem{T, U}}, n_diff::Integer) where {T, U}
    append!(
        vrss,
        fill(
            ValRootSystem(
                ValRoot(true, RootEdge{T, U}[], U[], U[]),
                ValRoot{T, U}[]
            ),
            n_diff
        )
    )

    return nothing
end

standard_segments(r::Root) = filter(!is_augmented, og_segment_id.(edges(r)))
standard_segments(vr::ValRoot) = og_segment_id.(edges(vr))

pixels(r::Root) = [(x, y) for (x, y) in zip(xs(r), ys(r))]
pixels(vr::ValRoot) = [(x, y) for (x, y) in zip(xs(vr), ys(vr))]