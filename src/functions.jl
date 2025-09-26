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
    create_facetsets(V ,face_boundary, cells, nodes)

Convert Comodo.jl face boundary to Ferrite.jl face boundary
"""
function create_facetsets(face_boundary, cells)
    strs = OrderedSet{FacetIndex}()
    
    for (element_id, element) in enumerate(cells)
        element_faces = Ferrite.faces(element)
        for (local_face_id, element_face) in enumerate(element_faces)
            # Check if this element face matches any boundary face
            for face_vertices in face_boundary
                if Set(element_face) == Set(face_vertices)
                    push!(strs, Ferrite.FacetIndex(element_id, local_face_id))
                    break  # Move to next element face
                end
            end
        end
    end
    
    return strs
end
###########################################################
###########################################################
function create_facetsets_twoD(face_boundary, cells)
    strs = OrderedSet{FacetIndex}()
    
    for (element_id, element) in enumerate(cells)
        element_faces = Ferrite.edges(element)
        for (local_face_id, element_face) in enumerate(element_faces)
            # Check if this element face matches any boundary face
            for face_vertices in face_boundary
                if Set(element_face) == Set(face_vertices)
                    push!(strs, Ferrite.FacetIndex(element_id, local_face_id))
                    break  # Move to next element face
                end
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
function ComodoToFerrite(F, V, ::Type{Ferrite.Quadrilateral}; Eb = nothing, Cb = nothing)
    cells = [Ferrite.Quadrilateral((e[1], e[2], e[3], e[4])) for e in F]
    nodes = [Ferrite.Node((e[1], e[2])) for e in V]

    if Eb === nothing && Cb === nothing
        return Grid(cells, nodes)
    else
        Fb_left = Eb[Cb.==4]
        Fb_bottom = Eb[Cb.==1]
        Fb_top = Eb[Cb.==3]
        Fb_right = Eb[Cb.==2]
        left = create_facetsets_twoD(Fb_left, cells)
        bottom = create_facetsets_twoD(Fb_bottom, cells)
        right = create_facetsets_twoD(Fb_right, cells)
        top = create_facetsets_twoD(Fb_top, cells)
        facetsets = Dict("left" => left, "bottom" => bottom, "right" => right, "top" => top)
        return Grid(cells, nodes, facetsets=facetsets)
    end 
    
end
###########################################################
###########################################################
"""
      ComodoToFerrite(F1, V1,Ferrite.Triangle)

Convert Tri3 of Comodo.jl to Ferrite.Triangle mesh in Ferrite.jl
"""
function ComodoToFerrite(F1, V1, ::Type{Ferrite.Triangle}, )

    cells = [Ferrite.Triangle((e[1], e[2], e[3])) for e in F1]
    nodes = [Ferrite.Node((e[1], e[2])) for e in V1]

    return Grid(cells, nodes)
end
###########################################################
###########################################################
function get_boundary_points(grid, facets, ::Type{Faces}, ::Type{T}) where {
    T<:Union{Ferrite.Quadrilateral,Ferrite.Triangle}
}
    facet_points = Point2f[]
    for facet in facets
        cell = grid.cells[facet[1]]
        facet_nodes = Ferrite.facets(cell)[facet[2]]
        for n in facet_nodes
            push!(facet_points, Point2f(Ferrite.get_node_coordinate(grid.nodes[n]).data...))
        end
    end
    facet_points = unique(facet_points)
    return facet_points
end
###########################################################
###########################################################

"""
     ComodoToFerrite(E, V, Ferrite.Hexahedron )

Convert Hex8 of Comodo.jl to Ferrite.Hexahedron mesh in Ferrite.jl
"""
function ComodoToFerrite(E, V, ::Type{Ferrite.Hexahedron}; Fb = nothing,  Cb= nothing)

    cells = [Ferrite.Hexahedron((e[1], e[2], e[3], e[4], e[5], e[6], e[7], e[8])) for e in E]
    nodes = [Ferrite.Node((e[1], e[2], e[3])) for e in V]

    if Fb === nothing && CFb_type === nothing
        return Grid(cells, nodes)
    else
        # based on Ferrite.jl
        Fb_bottom = Fb[Cb.==1]  # Bottom face (1)
        Fb_front = Fb[Cb.==3]   # Front face (2)
        Fb_top = Fb[Cb.==2]     # Top face (6)
        Fb_back = Fb[Cb.==4]    # Back face (4)
        Fb_right = Fb[Cb.==5]   # Right face (3)
        Fb_left = Fb[Cb.==6]    # Left face (5)

        left = create_facetsets( Fb_left, cells)
        bottom = create_facetsets(Fb_bottom, cells)
        right = create_facetsets( Fb_right, cells)
        back = create_facetsets(Fb_back, cells,)
        top = create_facetsets(Fb_top, cells,)
        front = create_facetsets(Fb_front, cells,)
        facetsets = Dict("left" => left, "bottom" => bottom, "right" => right, "back" => back, "top" => top, "front" => front)
        return Grid(cells, nodes, facetsets=facetsets)
    end
end

###########################################################
###########################################################
"""
      ComodoToFerrite(E, V, Ferrite.Tetrahedron)

Convert Tet4 of Comodo.jl to Ferrite.Tetrahedron mesh in Ferrite.jl
"""
function ComodoToFerrite(E, V, ::Type{Ferrite.Tetrahedron}; Fb = nothing,  Cb = nothing )

    cells = [Ferrite.Tetrahedron((e[1], e[2], e[3], e[4])) for e in E]
    nodes = [Ferrite.Node((e[1], e[2], e[3])) for e in V]

    if Fb === nothing && Cb === nothing
        return Grid(cells, nodes)
    else
        # based on Ferrite.jl
        Fb_bottom = Fb[Cb.==1]  # Bottom face (1)
        Fb_front = Fb[Cb.==3]   # Front face (2)
        Fb_top = Fb[Cb.==2]     # Top face (6)
        Fb_back = Fb[Cb.==4]    # Back face (4)
        Fb_right = Fb[Cb.==6]   # Right face (3)
        Fb_left = Fb[Cb.==5]    # Left face (5)

        left = create_facetsets(Fb_left, cells)
        bottom = create_facetsets(Fb_bottom, cells)
        right = create_facetsets(Fb_right, cells)
        back = create_facetsets(Fb_back, cells)
        top = create_facetsets(Fb_top, cells)
        front = create_facetsets(Fb_front, cells)
        facetsets = Dict("left" => left, "bottom" => bottom, "right" => right, "back" => back, "top" => top, "front" => front)
        return Grid(cells, nodes, facetsets=facetsets)
    end
    return 
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

    numElements = Ferrite.getncells(grid)
    CF_type = repeat(1:6,numElements) # Allocate face color/label data

    F = element2faces(E)

    F_uni, indUni ,c_uni = gunique(F,return_index=Val(true),return_counts=Val(true),sort_entries=true)
    Lb = isone.(c_uni)
    Fb = F_uni[Lb]
    CF_type_uni = CF_type[indUni]
    CFb_type = CF_type_uni[Lb]

    return E , V, F, Fb, CFb_type
end
###########################################################
###########################################################
function FerriteToComodo(grid, ::Type{Ferrite.Tetrahedron})
    E = Vector{Tet4{Int64}}()
    for i in eachindex(grid.cells)
        push!(E, grid.cells[i].nodes)
    end
    V = [Point{3,Float64}(Ferrite.get_node_coordinate(n)...) for n in grid.nodes]

    numElements = Ferrite.getncells(grid)
    CF_type = repeat(1:6,numElements) # Allocate face color/label data
    
    F = element2faces(E)

    F_uni, indUni, c_uni = gunique(F, return_index=Val(true), return_counts=Val(true), sort_entries=true)
    Lb = isone.(c_uni)
    Fb = F_uni[Lb]
    CF_type_uni = CF_type[indUni]
    CFb_type = CF_type_uni[Lb]
    return E, V, F, Fb, CFb_type
end