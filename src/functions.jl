"""
Functions to convert Comodo.jl Mesh to Ferrite.jl Mesh
In Ferrite.jl, Grid contains information about geometry, mesh, faces (edge 2D cases)
Supported Mesh type in Ferrite.jl:
Line
QuadraticLine
Quadrilateral 
QuadraticQuadrilateral
Hexahedron
Wedge
Pyramid
SerendipityQuadraticHexahedron
Triangle
QuadraticTriangle
Tetrahedron

In Ferrite.jl, function "generate_grid" generates Grid
It returns: Grid(cells, nodes, facetsets = facetsets)
So we need a function to get cells, nodes from Comodo.jl and then convert it 
to the type of Ferrite.jl. 
facetsets = facetsets: are face set in Ferrite

Name of fucntion:
CTF:  ComodoToFerrite (can be changed)
"""
###########################################################
###########################################################
"""
      ComodoToFerrite(F1, V1, Ferrite.Quadrilateral)

Convert Quad4 of Comodo.jl to Ferrite.Quadrilateral mesh in Ferrite.jl

"""
function ComodoToFerrite(F1, V1, ::Type{Ferrite.Quadrilateral})

    cells = [Ferrite.Quadrilateral((e[1], e[2], e[3], e[4])) for e in F1]
    nodes = [Ferrite.Node((e[1], e[2])) for e in V1]

    return Grid(cells, nodes)
end

###########################################################
###########################################################
"""
      ComodoToFerrite(F1, V1,Ferrite.Triangle)

Convert Tri3 of Comodo.jl to Ferrite.Triangle mesh in Ferrite.jl
"""
function ComodoToFerrite(F1, V1, ::Type{Ferrite.Triangle})

    cells = [Ferrite.Triangle((e[1], e[2], e[3])) for e in F1]
    nodes = [Ferrite.Node((e[1], e[2])) for e in V1]

    return Grid(cells, nodes)
end
###########################################################
###########################################################
"""
     ComodoToFerrite(E, V, Ferrite.Hexahedron )

Convert Hex8 of Comodo.jl to Ferrite.Hexahedron mesh in Ferrite.jl
"""
function ComodoToFerrite(E, V, ::Type{Ferrite.Hexahedron})

    cells = [Ferrite.Hexahedron((e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8])) for e in E]
    nodes = [Ferrite.Node((e[1], e[2], e[3])) for e in V]

    return Grid(cells, nodes)
end
###########################################################
###########################################################
"""
      ComodoToFerrite(E, V, Ferrite.Tetrahedron)

Convert Tet4 of Comodo.jl to Ferrite.Tetrahedron mesh in Ferrite.jl
"""
function ComodoToFerrite(E, V, ::Type{Ferrite.Tetrahedron})

    cells = [Ferrite.Tetrahedron((e[1], e[2], e[3], e[4])) for e in E]
    nodes = [Ferrite.Node((e[1], e[2], e[3])) for e in V]

    return Grid(cells, nodes)
end