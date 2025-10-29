module ComodoFerrite

using Comodo
using Ferrite
using OrderedCollections
using GeometryBasics
using GLMakie

struct Faces end
struct Nodes end

export Faces, Nodes
export Comodo

#################

### ComodoToFerrite
export create_facetsets
export create_facetsets_twoD
export ComodoToFerrite 
export get_boundary_points

### FerriteToComodo functions
export FerriteToComodo

include("functions.jl")

end
