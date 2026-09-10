module RootUntangling

export RootGraph, get_rootgraph # data to graph
export get_subgraphs # separate graphs
export solve_rsa # solving
export Root, get_roots, get_root_systems, examine, curve_length # get and examine roots
export greedy_switch, roughness # postprocess roots
export graphplot, rootplot, annotation_plot # plotting
export get_annotation_dict, write_annotation # write output
# debugging exports
export read_data, get_edge_info, get_vertex_info # data reading
export get_pregraph # data to graph
export src, dst, vertices, rootedge, segment_id, width, pred_primary, is_augmented, are_connected # edges
export id, edges, x, y, pred_split, xs, ys, rootvertex # vertices
export V₀, V₊, V, E₀, E₊, E, E₂, neighbors # graphs
export get_re_classification_dict, get_rv_classification_dict # get results

using JuMP, Statistics, GLMakie, GLMakie.Colors

include("data_reading.jl"); # read in data to dictionaries
include("data_cleaning.jl"); # clean dictionaries to expected format for vertex and edge info
include("pregraph_types.jl"); # data types for preliminary graph (segments from scan not yet divided into multiple possible roots)
include("pregraph_construction.jl"); # construct preliminary graph from data
include("rootgraph_types.jl"); # data types for rootgraph
include("rootgraph_functions.jl"); # rootgraph functions needed for the model
include("rootgraph_construction.jl"); # construct rootgraph from preliminary graph
include("rootgraph_clustering.jl"); # cluster disconnected graphs in a rootgraph
include("solving.jl"); # solve problem based on rootgraph
include("classification_extraction.jl"); # get classifications from model
include("root_types.jl"); # type to represent resulting roots
include("root_construction.jl"); # get the roots
include("root_postprocessing.jl"); # postprocessing on roots
include("plotting.jl"); # visualisation
include("annotation.jl"); # annotate graphs with root identities

end
