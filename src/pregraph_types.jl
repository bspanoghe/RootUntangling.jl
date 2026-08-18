# place holder datatypes for edges, vertices and graphs

isspecial(v::Integer) = v < 0

mutable struct Segment{T, U}
    id::T
    vertices::Vector{T}
    width::U
    pred_primary::Union{U, Missing} #! only for testing first 7 ROI - remove later
    xs::Vector{<:Number}
    ys::Vector{<:Number}
end
Segment(
    id::T, vertices::Vector{T}, width::U,
    pred_primary::Union{U, Missing}, xs, ys
) where {T, U} = (
    Segment(id, vertices, width, pred_primary, xs, ys)
)

Segment(id, vertices) = Segment(id, vertices, NaN, NaN, Number[], Number[])
id(s::Segment) = s.id
vertices(s::Segment) = s.vertices
width(s::Segment) = s.width
pred_primary(s::Segment) = s.pred_primary
xs(s::Segment) = s.xs
ys(s::Segment) = s.ys

isspecial(s::Segment) = any(isspecial(v) for v in vertices(s))

mutable struct MetaVertex{T, U}
    id::T
    name::Symbol
    x::U
    y::U
    pred_split::Union{U, Missing} #! only for testing first 7 ROI - remove later
end
function MetaVertex(id::T, name::Symbol, x::U, y::U) where {T, U}
    return MetaVertex(id, name, x, y, NaN)
end
function MetaVertex(id::T, name::Symbol) where {T}
    return MetaVertex(id, name, NaN, NaN)
end

id(mv::MetaVertex) = mv.id
name(mv::MetaVertex) = mv.name
x(mv::MetaVertex) = mv.x
y(mv::MetaVertex) = mv.y
pred_split(mv::MetaVertex) = mv.pred_split

isspecial(mv::MetaVertex) = id(mv) < 0
coords(mv::MetaVertex) = [x(mv), y(mv)]
distance(mv1::MetaVertex, mv2::MetaVertex) = (coords(mv1) - coords(mv2)) .^ 2 |> sum |> sqrt


struct PreGraph{T, U, V}
    vertices::Vector{T}
    segments::Vector{Segment{T, U}}
    metavertexdict::Dict{T, MetaVertex{T, V}}
    neighbordict::Dict{T, Vector{T}}

    # check that ids of normal vertices equals `1:length(normal_vertices)` and ids of special vertices equals `-1:-length(special_vertices)`
    function PreGraph(
            vertices::Vector{T}, segments::Vector{Segment{T, U}},
            metavertexdict::Dict{T, MetaVertex{T, V}}, neighbordict::Dict{T, Vector{T}}
        ) where {T, U, V}

        normalvertices = [v for v in vertices if !isspecial(v)]
        specialvertices = [v for v in vertices if isspecial(v)]
        @assert issetequal(normalvertices, 1:length(normalvertices))
        @assert issetequal(specialvertices, -1:-1:-length(specialvertices))

        return new{T, U, V}(vertices, segments, metavertexdict, neighbordict)
    end
end

# automatically generate neighbordict from other fields
function PreGraph(_vertices::Vector{T}, _segments::Vector{Segment{T, U}}, _metavertexdict::Dict{T, MetaVertex{T, V}}) where {T, U, V}
    neighbordict = Dict{T, Vector{T}}()
    for s in _segments
        v1, v2 = vertices(s)
        v1_nbs = get(neighbordict, v1, T[])
        v2_nbs = get(neighbordict, v2, T[])
        v2 in v1_nbs || (neighbordict[v1] = [v1_nbs; v2])
        v1 in v2_nbs || (neighbordict[v2] = [v2_nbs; v1])
    end
    return PreGraph(_vertices, _segments, _metavertexdict, neighbordict)
end

vertices(pg::PreGraph) = pg.vertices
segments(pg::PreGraph) = pg.segments
segments(pg::PreGraph{T, U}, v::T) where {T, U} = [s for s in segments(pg) if v in vertices(s)]
segments(pg::PreGraph{T, U}, mv::MetaVertex{T, U}) where {T, U} = segments(pg, id(mv))
getmetavertex(pg::PreGraph{T, U}, v::T) where {T, U} = pg.metavertexdict[v]
getmetavertices(pg::PreGraph) = [getmetavertex(pg, v) for v in vertices(pg)]
neighbors(pg::PreGraph{T, U, V}, v::T) where {T, U, V} = pg.neighbordict[v]

outer_vertices(pg::PreGraph) = [v for v in vertices(pg) if length(neighbors(pg, v)) == 1]
inner_vertices(pg::PreGraph) = [v for v in vertices(pg) if length(neighbors(pg, v)) > 1]
