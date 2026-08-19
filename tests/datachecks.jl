using Plots
using Pkg; Pkg.activate(".")
using Graphs

include(pwd() * "/src/graph_types.jl");
include(pwd() * "/src/data_reading.jl");
include(pwd() * "/src/data_cleaning.jl");
include(pwd() * "/src/plotting.jl");

cd(@__DIR__)

open("datachecks/data_output.txt", "w") do f
    write(f, "")
end

for roi_nr in 1:7
    println("Building my ROI $(roi_nr)!")

    begin
        edgelines = readlines("../clouddata/branchpoints/ROI_$(roi_nr)/segment_info_with_coords.csv")
        edge_data_dict = process_data(edgelines, :Segment_ID)

        vertexlines = readlines("../clouddata/branchpoints/ROI_$(roi_nr)/bp1_segments_grouped.csv")
        vertex_data_dict = process_data(vertexlines, :Node)

        # clean
        _branchpoints = getbranchpoints(vertex_data_dict)
        _segments = getsegments(_branchpoints, vertex_data_dict, edge_data_dict)

        single_vertex_segments = remove_single_vertex_segments!(_segments, debug_output = true)
        merged_branchpoints = cluster_vertices!(_segments, _branchpoints, dist_threshold = 3, debug_output = true)

        # turn into graph
        g = PreGraph(_branchpoints, _segments)
    end

    bad_vertices = [
        [getbranchpoint(g, bv) for bv in reduce(vcat, vertices.(single_vertex_segments), init = [])];
        merged_branchpoints
    ]
    open("datachecks/data_output.txt", "a") do f
        warning_lines = [
            isempty(merged_branchpoints) ? "" : "The following cluster vertices were made: $(name.(merged_branchpoints))";
            [
                "Edge $(id(s)) connected to vertices $(
                        [name(getbranchpoint(g, v)) * " ($(x(getbranchpoint(g, v)))/$(y(getbranchpoint(g, v))))" for v in vertices(s)]
                    ).\n"
                    for s in single_vertex_segments
            ]
        ]

        write(f, "ROI$(roi_nr)\n", warning_lines..., "\n")
    end

    # plot
    p = plot(
        label = false, yflip = true, aspect_ratio = :equal, size = (2000, 2000)
    )
    plot!(
        [v for v in getbranchpoints(g) if v.name in name.(bad_vertices)],
        markerstrokewidth = 0, markersize = 3, color = :fuchsia, label = false
    )
    plot!(
        [v for v in getbranchpoints(g) if !(v.name in name.(bad_vertices))],
        markerstrokewidth = 0, markersize = 3, color = :cyan, label = false
    )
    for e in edges(g)
        plot!(e, g, linecolor = :teal, label = false)
    end

    savefig("datachecks/roi$(roi_nr).svg")
end
