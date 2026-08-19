module RootUntangling

export RootGraph, get_rootgraph # data to graph
export get_subgraphs # separate graphs
export solve_rsa # solving
export Root, get_roots, examine, curve_length # get and examine roots
export greedy_switch, tortuosity # postprocess roots
export hypothesis_plot # fancy plotting
# debugging exports
export read_data, get_edge_info, get_vertex_info # data reading
export get_pregraph # data to graph
export src, dst, vertices, width, pred_primary, is_augmented, are_connected # edges
export id, edges, x, y, pred_split, xs, ys # vertices
export V₀, V₊, V, E₀, E₊, E, E₂, segments₀, segments₊, segments, neighbors # graphs
export get_edge_classification_dict, get_segment_classification_dict, get_vertex_classification_dict # get results

using JuMP, Statistics, Plots

include("data_reading.jl"); # read in data to dictionaries
include("data_cleaning.jl"); # clean dictionaries to expected format for vertex and edge info
include("pregraph_types.jl"); # data types for preliminary graph (segments from scan not yet divided into multiple possible roots)
include("pregraph_construction.jl"); # construct preliminary graph from data
include("rootgraph_types.jl"); # data types for root graph
include("rootgraph_functions.jl"); # root graph functions needed for the model
include("rootgraph_construction.jl"); # construct root graph from preliminary graph
include("rootgraph_clustering.jl"); # cluster disconnected graphs in a root graph
include("solving.jl"); # solve problem based on root graph
include("classification_extraction.jl") # get classifications from model
include("root_types.jl") # type to represent resulting roots
include("root_construction.jl") # get the roots
include("root_postprocessing.jl") # postprocessing on roots
include("plotting.jl"); # visualisation

end
