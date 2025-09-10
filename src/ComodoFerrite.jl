module ComodoFerrite
using Ferrite
using OrderedCollections
using GeometryBasics
struct Faces end
struct Nodes end

export Faces, Nodes

#################

export create_facetsets
export ComodoToFerrite 
export get_boundary_points
include("functions.jl")
end
