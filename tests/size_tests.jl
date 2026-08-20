using Pkg; Pkg.activate("./scripts")
using RootUntangling, RootUntangling.Plots

file = "tests/outputs/sizes.csv"
dist_threshold = 3
pₛ_values = [0.1, 0.5]

begin
    open(file, "w") do f
        write(f, "roi, p_s, n_v, n_e, n_c")
    end

    for pₛ in pₛ_values
        for roi_nr in 1:7
            filename_segments = "./data/ROI_$(roi_nr)/segment_info_with_coords.csv"
            filename_vertices = "./data/ROI_$(roi_nr)/bp1_segments_grouped.csv"

            rg = get_rootgraph(filename_segments, filename_vertices; dist_threshold, pₛ, reverse_y = true)

            connections_by_vertex = [
                [
                        [edges(v)[i], edges(v)[j]]
                        for i in eachindex(edges(v)) for j in eachindex(edges(v))
                        if i > j
                    ]
                    for v in V₀(rg)
            ]
            connections = reduce(vcat, connections_by_vertex)

            n_v = length(V(rg))
            n_e = length(E(rg))
            n_c = length(connections)

            open(file, "a") do f
                write(f, "\n$roi_nr, $pₛ, $n_v, $n_e, $n_c")
            end
        end
    end
end

datadict = (
    readlines(file) .|>
        (x -> split(x, ", ")) |>
        (x -> reduce(hcat, x)) |>
        permutedims |>
        x -> [Symbol(c[1]) => parse.(Float64, c[2:end]) for c in eachcol(x)] |>
        Dict
)

begin
    p_vc2 = scatter(
        datadict[:n_v][datadict[:p_s] .== pₛ_values[1]], datadict[:n_c][datadict[:p_s] .== pₛ_values[1]],
        legend = false, title = "pₛ = $(pₛ_values[1])"
    )
    p_vc3 = scatter(
        datadict[:n_v][datadict[:p_s] .== pₛ_values[2]], datadict[:n_c][datadict[:p_s] .== pₛ_values[2]],
        legend = false, title = "pₛ = $(pₛ_values[2])", xlabel = "Number of vertices"
    )
    plot(p_vc2, p_vc3, layout = (2, 1), ylabel = "Number of connections")
end

begin
    p_ve2 = scatter(
        datadict[:n_v][datadict[:p_s] .== pₛ_values[1]], datadict[:n_e][datadict[:p_s] .== pₛ_values[1]],
        legend = false, title = "pₛ = $(pₛ_values[1])"
    )
    p_ve3 = scatter(
        datadict[:n_v][datadict[:p_s] .== pₛ_values[2]], datadict[:n_e][datadict[:p_s] .== pₛ_values[2]],
        legend = false, title = "pₛ = $(pₛ_values[2])", xlabel = "Number of vertices"
    )
    plot(p_ve2, p_ve3, layout = (2, 1), ylabel = "Number of edges")
end
