module RootUntangling

export SuperGraph, get_supergraph # data to graph
export get_subgraphs # separate graphs
export solve_rsa # solving
export Root, get_roots, examine, curve_length # get and examine roots
export greedy_switch, tortuosity # postprocess roots
export hypothesis_plot # fancy plotting
# debugging exports
export read_data, get_edge_info, get_vertex_info # data reading
export get_pregraph # data to graph
export src, dst, vertices, hyperedge, width, pred_primary, is_augmented, are_connected # edges
export id, edges, x, y, pred_split, xs, ys, hypervertex, exclusion_sets # vertices
export Vₕ₀, Vₕ₊, V₀, V₊, Vₕ, V, Eₕ₀, Eₕ₊, E₀, E₊, Eₕ, E, E₂, neighbors # graphs
export get_he_classification_dict, get_se_classification_dict, get_hv_classification_dict, get_sv_classification_dict # get results

using JuMP, Statistics, Plots

include("data_reading.jl"); # read in data to dictionaries
include("data_cleaning.jl"); # clean dictionaries to expected format for vertex and edge info
include("pregraph_types.jl"); # data types for preliminary graph (segments from scan not yet divided into multiple possible roots)
include("pregraph_construction.jl"); # construct preliminary graph from data
include("supergraph_types.jl"); # data types for supergraph
include("supergraph_functions.jl"); # supergraph functions needed for the model
include("supergraph_construction.jl"); # construct supergraph from preliminary graph
include("supergraph_clustering.jl"); # cluster disconnected graphs in a supergraph
include("solving.jl"); # solve problem based on supergraph
include("classification_extraction.jl") # get classifications from model
include("root_types.jl") # type to represent resulting roots
include("root_construction.jl") # get the roots
include("root_postprocessing.jl") # postprocessing on roots
include("plotting.jl"); # visualisation

end
