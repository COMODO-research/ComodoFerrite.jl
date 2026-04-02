module ComodoFerrite

using Comodo
using Ferrite
using OrderedCollections
using GeometryBasics
using GLMakie

export Comodo
export Ferrite
#################

export addface!
export ComodoToFerrite 
export get_boundary_points
export faceset_to_cellid_faceid 
export FerriteToComodo
### FerriteToComodo functions
export FerriteToComodo
export boundary_facets
include("functions.jl")

end
