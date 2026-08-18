# get required information for vertices and segments

# ## vertex information
function get_vertex_info(datadict::Dict; segment_ids_colname::Symbol, x_colname::Symbol,
    y_colname::Symbol, lateral_score_colname::Symbol, y_transform::Function)

    ymax = [vertex_data[y_colname] for vertex_data in values(datadict)] |> maximum

    vertex_data_dict = [
        Symbol(entry[1]) => Dict([
            :segment_ids => entry[2][segment_ids_colname],
            :x => entry[2][x_colname],
            :y => y_transform(entry[2][y_colname]),
            :pred_split => get(entry[2], lateral_score_colname, missing) #! only for testing first 7 ROI - remove default `missing` later
        ])
        for entry in datadict
    ] |> Dict

    return vertex_data_dict
end

# ## segment information

function get_edge_info(datadict::Dict; dist_colname::Symbol,
    primary_score_colname::Symbol, coords_colname::Symbol, y_transform::Function)

    edge_data_dict = [
        entry[1] => Dict([
            :width => entry[2][dist_colname],
            :pred_primary => (
                get(entry[2], primary_score_colname, missing) |> x -> isequal(x, "") ? missing : x
            ), #! only for testing first 7 ROI - remove default `missing` later
            :xs => entry[2][coords_colname] .|> (x -> split(x, '/')) .|> last .|> x -> parse(Int64, x),
            :ys => entry[2][coords_colname] .|> (x -> split(x, '/')) .|> first .|> (x -> parse(Int64, x)) .|> y_transform
        ])
        for entry in datadict
    ] |> Dict

    return edge_data_dict
end