module ComodoFerrite
using Comodo
using Ferrite
using OrderedCollections
using GeometryBasics
struct Faces end
struct Nodes end

export Faces, Nodes

#################

### ComodToFerrite
export create_facetsets
export ComodoToFerrite 
export get_boundary_points


### FerriteToComodo functions
export FerriteToComodo
include("functions.jl")
end
