# read in data from file

function read_data(
        file::String, id_colname::Symbol;
        delim::Char = ',', groupdelim::Char = '"', superdelim::Char = ';'
    )

    lines = readlines(file)

    processed_lines = lines .|>
        x -> split_line(x, delim, groupdelim) .|>
        x -> autoparse(x, delim, superdelim)

    header = processed_lines[1] .|> Symbol
    data = processed_lines[2:end] |> x -> reduce(hcat, x)

    data_dict = Dict(
        [header[i] => vectorpromote(data[i, :]) for i in eachindex(header)]
    )

    idwise_dict = Dict(
        [
            id => Dict([colname => data_dict[colname][i] for colname in header if colname != id_colname])
                for (i, id) in enumerate(data_dict[id_colname])
        ]
    )

    return idwise_dict
end

function split_line(line, delim::Char, groupdelim::Char)
    line_groups = line |>
        x -> split(x, groupdelim) |>
        x -> x[.!isempty.(x)] |>
        x -> single_strip.(x, delim)

    line_elements = [
        line_group |>
            x -> isodd(i) ? split(x, delim) : [x]
            for (i, line_group) in enumerate(line_groups)
    ] |> x -> reduce(vcat, x)

    return line_elements
end

"""
    single_strip(s::AbstractString, chars::Vector{Char})

`strip` but with a maximum recursive depth of 1.
"""
single_strip(s::AbstractString, chars::Vector{Char}) = [
    c for (i, c) in enumerate(s) if !(c in chars && (i == 1 || i == length(s)))
] |> x -> *(x...)
single_strip(s::AbstractString, char::Char) = single_strip(s, [char])

function autoparse(line::AbstractString, delim::Char, superdelim::Char)
    # try parsing as number or bool
    !isnothing(tryparse(Int64, line)) && return parse(Int64, line)
    !isnothing(tryparse(Float64, line)) && return parse(Float64, line)
    !isnothing(tryparse(Bool, lowercase(line))) && return parse(Bool, lowercase(line))

    # don't change anything if no delimiters are present
    !(delim in line) && !(superdelim in line) && return line

    # line with delimiters gets split into elements, each of which are attempted to be parsed again
    line_elements = split(line, superdelim) |>
        x -> split.(x, delim) .|>
        x -> strip.(x, [[' ', '(', ')']])

    return [unwrap([autoparse(_le, delim, superdelim) for _le in le]) for le in line_elements] |> unwrap
end

unwrap(x::AbstractVector) = length(x) == 1 ? unwrap(x[1]) : x
unwrap(x) = x

function vectorpromote(xs)
    if any(isa.(xs, Vector))
        xs = wrap.(xs)
    end
    xs = collect(promote(xs...))
    return xs
end

wrap(x) = [x]
wrap(x::Vector) = x
