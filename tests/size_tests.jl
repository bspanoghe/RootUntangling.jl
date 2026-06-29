using Pkg; Pkg.activate("./scripts")
using RootUntangling, RootUntangling.Plots

file = "tests/outputs/sizes.csv"
dist_threshold = 3

begin
    open(file, "w") do f
        write(f, "roi, n_h, n_v, n_e, n_c")
    end

    for n_h in 2:3
        for roi_nr in 1:7
            filename_segments = "./clouddata/branchpoints/ROI_$(roi_nr)/segment_info_with_coords.csv";
            filename_vertices = "./clouddata/branchpoints/ROI_$(roi_nr)/bp1_segments_grouped.csv";

            sg = get_supergraph(filename_segments, filename_vertices; dist_threshold, num_hypotheses = n_h)

            connections_by_vertex = [
                [
                    [edges(v)[i], edges(v)[j]] 
                    for i in eachindex(edges(v)) for j in eachindex(edges(v))
                    if i > j
                ]
                for v in V₀(sg)
            ]
            connections = reduce(vcat, connections_by_vertex)

            n_v = length(V(sg))
            n_e = length(E(sg))
            n_c = length(connections)

            open(file, "a") do f
                write(f, "\n$roi_nr, $n_h, $n_v, $n_e, $n_c")
            end
        end
    end
end

datadict = (
    readlines(file) .|> 
    (x -> split(x, ", ")) |>
    (x -> reduce(hcat, x)) |>
    permutedims |>
    x -> [Symbol(c[1]) => parse.(Int64, c[2:end]) for c in eachcol(x)] |>
    Dict
)

begin
    p_vc2 = scatter(datadict[:n_v][datadict[:n_h] .== 2], datadict[:n_c][datadict[:n_h] .== 2],
        legend = false, title = "2 hypotheses"
    )
    p_vc3 = scatter(datadict[:n_v][datadict[:n_h] .== 3], datadict[:n_c][datadict[:n_h] .== 3],
        legend = false, title = "3 hypotheses", xlabel = "Number of vertices"
    )
    plot(p_vc2, p_vc3, layout = (2, 1), ylabel = "Number of connections")
end

begin
    p_ve2 = scatter(datadict[:n_v][datadict[:n_h] .== 2], datadict[:n_e][datadict[:n_h] .== 2],
        legend = false, title = "2 hypotheses"
    )
    p_ve3 = scatter(datadict[:n_v][datadict[:n_h] .== 3], datadict[:n_e][datadict[:n_h] .== 3],
        legend = false, title = "3 hypotheses", xlabel = "Number of vertices"
    )
    plot(p_ve2, p_ve3, layout = (2, 1), ylabel = "Number of edges")
end