# # edges
abstract type AbstractEdge{T} end
src(ae::AbstractEdge) = ae.src
dst(ae::AbstractEdge) = ae.dst
vertices(ae::AbstractEdge) = (src(ae), dst(ae))

is_augmented(ae::AbstractEdge) = any(is_augmented.(vertices(ae)))
Base.:(==)(ae1::AbstractEdge, ae2::AbstractEdge) = vertices(ae1) == vertices(ae2)
Base.unique(aes::Vector{<:AbstractEdge}) = unique(x -> vertices(x), aes) # doesn't automatically use my equality operator :(
Base.unique!(aes::Vector{<:AbstractEdge}) = unique!(x -> vertices(x), aes) # sad times

Base.show(io::IO, ae::AbstractEdge) = print(io, "$(typeof(ae).name.name)$(vertices(ae))")
Base.show(io::IO, aes::Vector{<:AbstractEdge}) = print(io, "$(typeof(aes).name.name)$(vertices.(aes))")

"""
    HyperEdge{T, U}

Represents a segment from a root scan, which may contain one or multiple roots.
"""
struct HyperEdge{T, U} <: AbstractEdge{T}
    src::T
    dst::T
    width::U
    pred_primary::Union{U, Missing} #! only for testing first 7 ROI - remove Missing Union later

    HyperEdge(src::T, dst::T, width::U, pred_primary::Union{U, Missing}) where {T, U} = (
        new{T, U}(sort([src, dst])..., width, pred_primary)
    )
end
width(he::HyperEdge) = he.width
pred_primary(he::HyperEdge) = he.pred_primary

HyperEdge(src::T, dst::T) where {T} = HyperEdge(src, dst, NaN, NaN)

"""
    SingularEdge{T, U}

Represents an edge in between two vertices of different hypervertices.
"""
struct SingularEdge{T, U} <: AbstractEdge{T}
    src::T
    dst::T
    hyperedge::HyperEdge{T, U}
    
    SingularEdge(src::T, dst::T, hyperedge::HyperEdge{T, U}) where {T, U} = (
        new{T, U}(sort([src, dst])..., hyperedge)
    )    
end
hyperedge(se::SingularEdge) = se.hyperedge

# vertices
abstract type AbstractVertex{T, U} end
id(av::AbstractVertex) = av.id
edges(av::AbstractVertex) = av.edges
edges(avs::Vector{<:AbstractVertex}) = reduce(vcat, edges.(avs), init = eltype(avs)[])

Base.show(io::IO, av::AbstractVertex) = print(io, "$(typeof(av).name.name)$( (id(av), edges(av)) )")
Base.show(io::IO, avs::Vector{<:AbstractVertex}) = print(io, "$(typeof(avs).name.name)$( id.(avs) )")

"""
    HyperVertex{T, U}

Represents a branchpoint or endpoint of a segment in the image.
"""
struct HyperVertex{T, U} <: AbstractVertex{T, U}
    id::T
    edges::Vector{HyperEdge{T, U}}
    x::U
    y::U
    pred_split::Union{U, Missing} #! only for testing first 7 ROI - remove later
    vertices::Vector{T}

    HyperVertex(id::T, edges::Vector{HyperEdge{T, U}}, x::U, y::U, pred_split::Union{U, Missing}, vertices::Vector{T}) where {T, U} = (
        new{T, U}(id, sort(edges, by = e -> src(e)), x, y, pred_split, sort(vertices))
    )
end

x(hv::HyperVertex) = hv.x
y(hv::HyperVertex) = hv.y
pred_split(hv::HyperVertex) = hv.pred_split
vertices(hv::HyperVertex) = hv.vertices

gethypervertex(v::T, Vₕ::Vector{HyperVertex{T, U}}) where {T, U} = Vₕ[findfirst(hv -> id(hv) == v, Vₕ)]
neighbor(hv::HyperVertex{T, U}, he::HyperEdge{T, U}, Vₕ::Vector{HyperVertex{T, U}}) where {T, U} = (
    vertices(he)[findfirst(v -> v != id(hv), vertices(he))] |> (v -> gethypervertex(v, Vₕ))
)
is_augmented(hv::HyperVertex) = any(is_augmented.(vertices(hv)))

coords(hv::HyperVertex) = (x(hv), y(hv))


"""
    SingularVertex{T, U}

Represents a point on a single root.
"""
struct SingularVertex{T, U} <: AbstractVertex{T, U}
    id::T
    edges::Vector{SingularEdge{T, U}}
    hypervertex::HyperVertex{T, U}

    SingularVertex(id::T, edges::Vector{SingularEdge{T, U}}, hypervertex::HyperVertex{T, U}) where {T, U} = (
        new{T, U}(id, sort(edges, by = e -> src(e)), hypervertex)
    )
end
hypervertex(sv::SingularVertex) = sv.hypervertex

is_augmented(v::T) where {T} = v < 0
is_augmented(sv::SingularVertex) = is_augmented(id(sv))
x(sv::SingularVertex) = x(hypervertex(sv))
y(sv::SingularVertex) = y(hypervertex(sv))
coords(sv::SingularVertex) = coords(hypervertex(sv))

# graphs
"""
    SuperGraph{T, U}

Represents one or more root systems.

The (singular) edges E represent the possible presence of a single root, 
while hyperedges Eₕ are groups of edges at the same position, representing a segment in the scan which may contain multiple roots.
"""
struct SuperGraph{T, U}
    Vₕ₀::Vector{HyperVertex{T, U}}
    Vₕ₊::Vector{HyperVertex{T, U}}
    V₀::Vector{SingularVertex{T, U}}
    V₊::Vector{SingularVertex{T, U}}

    Eₕ₀::Vector{HyperEdge{T}}
    Eₕ₊::Vector{HyperEdge{T}}
    E₀::Vector{SingularEdge{T, U}} 
    E₊::Vector{SingularEdge{T, U}}
    function SuperGraph(Vₕ₀::Vector{HyperVertex{T, U}}, Vₕ₊::Vector{HyperVertex{T, U}},
        V₀::Vector{SingularVertex{T, U}}, V₊::Vector{SingularVertex{T, U}}) where {T, U}
        
        # sort all vertices
        sort!(Vₕ₀, by = x -> id(x))
        sort!(Vₕ₊, by = x -> -id(x)) # special vertices use negative integers as id
        sort!(V₀, by = x -> id(x))
        sort!(V₊, by = x -> -id(x)) # special vertices use negative integers as id

        # assert all ids are consecutive integers starting from 1 / -1
        @assert id.(Vₕ₀) == 1:length(Vₕ₀)
        @assert id.(Vₕ₊) == -1:-1:-length(Vₕ₊)
        @assert id.(V₀) == 1:length(V₀)
        @assert id.(V₊) == -1:-1:-length(V₊)

        # get edges and sort
        Eₕ₀ = [e for e in edges(Vₕ₀) if (!is_augmented(e))] |> unique |> es -> sort(es, by = e -> src(e))
        Eₕ₊ = edges(Vₕ₊) |> unique |> es -> sort(es, by = e -> src(e))
        E₀ = [e for e in edges(V₀) if !(is_augmented(e))] |> unique |> es -> sort(es, by = e -> src(e))
        E₊ = edges(V₊) |> unique |> es -> sort(es, by = e -> src(e))

        return new{T, U}(Vₕ₀, Vₕ₊, V₀, V₊, Eₕ₀, Eₕ₊, E₀, E₊)
    end
end
Vₕ₀(sg::SuperGraph) = sg.Vₕ₀
Vₕ₊(sg::SuperGraph) = sg.Vₕ₊
V₀(sg::SuperGraph) = sg.V₀
V₊(sg::SuperGraph) = sg.V₊

Eₕ₀(sg::SuperGraph) = sg.Eₕ₀
Eₕ₊(sg::SuperGraph) = sg.Eₕ₊
E₀(sg::SuperGraph) = sg.E₀
E₊(sg::SuperGraph) = sg.E₊

Vₕ(sg::SuperGraph) = [Vₕ₊(sg); Vₕ₀(sg)]
V(sg::SuperGraph) = [V₊(sg); V₀(sg)]
Eₕ(sg::SuperGraph) = [Eₕ₊(sg); Eₕ₀(sg)]
E(sg::SuperGraph) = [E₊(sg); E₀(sg)]

E₂(sv::SingularVertex) = [
    [edges(sv)[i], edges(sv)[j]] 
    for i in eachindex(edges(sv)) for j in eachindex(edges(sv))
    if (i > j) && !all(is_augmented.(edges(sv)[[i, j]]))
]
E₂(sg::SuperGraph) = E₂.(V₀(sg)) |> x -> reduce(vcat, x, init = eltype(x)[])
E₂(sv::SingularVertex, se::SingularEdge) = [c for c in E₂(sv) if se in c]

Base.length(sg::SuperGraph) = length(Vₕ₀(sg))

getsingularvertex(sg::SuperGraph{T, U}, v::T) where {T, U} = v > 0 ? V₀(sg)[v] : V₊(sg)[-v]
gethypervertex(sg::SuperGraph{T, U}, v::T) where {T, U} = v > 0 ? Vₕ₀(sg)[v] : Vₕ₊(sg)[-v]
width(sg::SuperGraph, se::SingularEdge) = width(hyperedge(se)) / minimum(order.(getsingularvertex.([sg], vertices(se))))

neighbor(sg::SuperGraph{T, U}, sv::SingularVertex{T, U}, se::SingularEdge{T}) where {T, U} = (
    vertices(se)[findfirst(v -> v != id(sv), vertices(se))] |> v -> getsingularvertex(sg, v)
)
neighbor(sg::SuperGraph{T, U}, hv::HyperVertex{T, U}, he::HyperEdge{T}) where {T, U} = (
    vertices(he)[findfirst(v -> v != id(hv), vertices(he))] |> v -> gethypervertex(sg, v)
)
neighbors(sg::SuperGraph{T, U}, av::AbstractVertex{T, U}) where {T, U} = [neighbor(sg, av, ae) for ae in edges(av)]

inner_vertices(sg::SuperGraph) = [v for v in V₀(sg) if length([n for n in neighbors(sg, Vₕ(v)) if !is_augmented(n)]) > 1]
outer_vertices(sg::SuperGraph) = [v for v in V₀(sg) if length([n for n in neighbors(sg, Vₕ(v)) if !is_augmented(n)]) == 1]

xs(sg::SuperGraph{T, U}, se::SingularEdge{T, U}) where {T, U} = x.(V(sg, se))
ys(sg::SuperGraph{T, U}, se::SingularEdge{T, U}) where {T, U} = y.(V(sg, se))

# some additional methods to make model definition look nicer
V(sg::SuperGraph{T, U}, se::SingularEdge{T, U}) where {T, U} = [getsingularvertex(sg, v) for v in vertices(se)]
V(sg::SuperGraph{T, U}, he::HyperEdge{T, U}) where {T, U} = [vertices(gethypervertex(v, Vₕ(sg))) for v in vertices(he)]
V(sg::SuperGraph{T, U}, hv::HyperVertex{T, U}) where {T, U} = [getsingularvertex(sg, v) for v in vertices(hv)]
Vₕ(sv::SingularVertex) = hypervertex(sv)
Vₕ(sg::SuperGraph, he::HyperEdge) = [gethypervertex(v, Vₕ(sg)) for v in [src(he), dst(he)]]

E(av::AbstractVertex) = edges(av)
E(avs::Vector{<:AbstractVertex}) = reduce(vcat, E.(avs), init = eltype(avs)[])
E(sg::SuperGraph{T, U}, hv::HyperVertex{T, U}) where {T, U} = E(V(sg, hv))
E(sg::SuperGraph{T, U}, he::HyperEdge{T, U}) where {T, U} = [e for e in E(sg) if hyperedge(e) == he] |> unique
E₀(av::AbstractVertex) = [e for e in E(av) if !is_augmented(e)]
Eₕ(hv::HyperVertex) = edges(hv)
Eₕ(sv::SingularVertex) = Eₕ(Vₕ(sv))
Eₕ(se::SingularEdge) = hyperedge(se)
Eₕ₀(sv::SingularVertex) = [he for he in Eₕ(sv) if !is_augmented(he)]