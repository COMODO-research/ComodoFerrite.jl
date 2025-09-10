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
     create_facetsets(E ,face_boundary, face_index::Int64)

Convert Comodo.jl face boundary to Ferrite.jl face boundary
"""
function create_facetsets(E ,face_boundary, face_index::Int64)
    strs = OrderedSet{FacetIndex}()
    for i in eachindex(Vector.(face_boundary))
        for j in eachindex(E)
            check = issubset(Set(face_boundary[i]), Set(E[j]))
            if check == true
               cell_index = j
               push!(strs, Ferrite.FacetIndex(cell_index,face_index) )
            end
        end
    end
    return strs
end
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
function ComodoToFerrite(E, V, Fb, CFb_type,  ::Type{Ferrite.Hexahedron})

    cells = [Ferrite.Hexahedron((e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8])) for e in E]
    nodes = [Ferrite.Node((e[1], e[2], e[3])) for e in V]
    # based on Ferrite.jl
    Fb_bottom = Fb[CFb_type.==1]  # Bottom face (1)
    Fb_front = Fb[CFb_type.==3]   # Front face (2)
    Fb_top = Fb[CFb_type.==2]     # Top face (6)
    Fb_back = Fb[CFb_type.==4]    # Back face (4)
    Fb_right = Fb[CFb_type.==5]   # Right face (3)
    Fb_left = Fb[CFb_type.==6]    # Left face (5)

    left = create_facetsets(E, Fb_left, 5)
    bottom = create_facetsets(E, Fb_bottom, 1)
    right = create_facetsets(E, Fb_right, 3)
    back = create_facetsets(E, Fb_back, 4)
    top = create_facetsets(E, Fb_top, 6)
    front = create_facetsets(E, Fb_front, 2)
    facetsets = Dict("left" => left, "bottom" => bottom, "right" => right, "back" => back, "top" => top, "front" => front)
    
    return Grid(cells, nodes, facetsets = facetsets)
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
###########################################################
###########################################################
"""
    get_boundary_points(grid, facets, Faces, Ferrite.Hexahedron)
    or 
    get_boundary_points(grid, facets, Faces, Ferrite.Tetrahedron)

function to convert FaceIndex boundary of Ferrite.jl to points for plotting the boundry condition
"""
function get_boundary_points(grid, facets, ::Type{Faces}, ::Type{T}) where {
    T<:Union{Ferrite.Hexahedron,Ferrite.Tetrahedron}
}
    facet_points = Point3f[]
    for facet in facets
        cell = grid.cells[facet[1]]
        facet_nodes = Ferrite.facets(cell)[facet[2]]
        for n in facet_nodes
            push!(facet_points, Point3f(Ferrite.get_node_coordinate(grid.nodes[n]).data...))
        end
    end
    facet_points = unique(facet_points)
    return facet_points
end
###########################################################
###########################################################
"""
     get_boundary_points(grid, facets, Faces, Ferrite.Hexahedron)
     or 
     get_boundary_points(grid, facets, Faces, Ferrite.Tetrahedron)

Convert a set of node indices from a Ferrite grid to 3D point coordinates for plotting.
"""
function get_boundary_points(grid, nodeset, ::Type{Nodes}, ::Type{T}) where {
    T<:Union{Ferrite.Hexahedron, Ferrite.Tetrahedron}
}
    nodesset = [Point{3, Float64}(Ferrite.get_node_coordinate(grid.nodes[i]).data...) for i in nodeset]
    return nodesset
end
###########################################################
###########################################################
"""

"""
function FerriteToComodo(grid, ::Type{Ferrite.Hexahedron})
     E = Vector{Hex8{Int64}}()
    for i in eachindex(grid.cells)
        push!(E , grid.cells[i].nodes)
    end
    V = [Point{3,Float64}(Ferrite.get_node_coordinate(n)...) for n in grid.nodes]

    F = element2faces(E)

    F_uni,_,c_uni = gunique(F,return_index=Val(true),return_counts=Val(true),sort_entries=true)
    Lb = isone.(c_uni)
    Fb = F_uni[Lb]
    return E , V, F, Fb
end