function addface!(grid, name::String, boundary_faces)
    # Create dictionaries for all possible face sizes
    face_dicts = Dict{Int, Dict{NTuple, Int}}()
    
    # Populate dictionaries, grouping by face size
    for (i, face) in enumerate(boundary_faces)
        sorted_face = Tuple(sort(collect(face)))
        n_vertices = length(sorted_face)
        
        # Create dictionary for this size if it doesn't exist
        if !haskey(face_dicts, n_vertices)
            face_dicts[n_vertices] = Dict{NTuple{n_vertices, Int}, Int}()
        end
        
        face_dicts[n_vertices][sorted_face] = i
    end
    
    # Get unique node indices from boundary faces
    boundary_nodes = Set{Int}()
    for face in boundary_faces
        union!(boundary_nodes, face)
    end
    
    # Find elements that contain any boundary node
    candidate_elements = Set{Int}()
    for (element_id, cell) in enumerate(grid.cells)
        cell_nodes = Ferrite.get_node_ids(cell)
        if !isdisjoint(cell_nodes, boundary_nodes)
            push!(candidate_elements, element_id)
        end
    end
    
    # Create new facetset - only loop through candidate elements
    new_facetset = Set{FacetIndex}()
    for element_id in candidate_elements
        cell = grid.cells[element_id]
        element_facets = Ferrite.facets(cell)
        for (local_facet_id, element_facet) in enumerate(element_facets)
            sorted_facet = Tuple(sort(collect(element_facet)))
            n_vertices = length(sorted_facet)
            
            # Check in the appropriate dictionary for this face size
            if haskey(face_dicts, n_vertices) && haskey(face_dicts[n_vertices], sorted_facet)
                push!(new_facetset, FacetIndex(element_id, local_facet_id))
            end
        end
    end
    
    # Add to grid's facetsets
    grid.facetsets[name] = new_facetset
    
    return grid
end

"""
      ComodoToFerrite(connectivty, V)

Convert Comodo.jl  mesh to Ferrite grid in Ferrite.jl
connectivty: for 3D is E and for 2D is F

"""
function ComodoToFerrite(connectivity, V)
    CellType = eltype(connectivity)
    
    FerriteCell = if CellType == Hex8{Int64}
        Ferrite.Hexahedron
    elseif CellType == Tet4{Int64}
        Ferrite.Tetrahedron
    elseif CellType == QuadFace{Int64}
        Ferrite.Quadrilateral
    elseif CellType == TriangleFace{Int64}
        Ferrite.Triangle
    else
        error("Unsupported cell type: $(CellType)")
    end
    
    # Determine dimension based on cell type
    dim = CellType in (QuadFace{Int64}, TriangleFace{Int64}) ? 2 : 3
    
    cells = [FerriteCell(Tuple(e)) for e in connectivity]
    nodes = [Ferrite.Node(ntuple(i -> v[i], dim)) for v in V]
    
    return Grid(cells, nodes)
end

function get_boundary_points(grid, dataset)
    # Determine dimension from grid nodes
    NodeType = eltype(grid.nodes)
    dim = NodeType.parameters[1]
    PointType = dim == 2 ? Point2f : Point3f
    
    # Extract points based on dataset type
    if dataset isa OrderedCollections.OrderedSet{Int64}
        # Direct node indices
        points = [PointType(Ferrite.get_node_coordinate(grid.nodes[i]).data...) 
                  for i in dataset]
    
    elseif dataset isa OrderedCollections.OrderedSet{FacetIndex}
        # Extract nodes from facets
        points = PointType[]
        for facet in dataset
            cell = grid.cells[facet[1]]
            facet_nodes = Ferrite.facets(cell)[facet[2]]
            for node_idx in facet_nodes
                push!(points, PointType(Ferrite.get_node_coordinate(grid.nodes[node_idx]).data...))
            end
        end
        unique!(points)  # In-place unique
    
    else
        error("Invalid dataset type: expected OrderedSet{Int64} or OrderedSet{FacetIndex}, got $(typeof(dataset))")
    end
    
    return points
end


"""

"""
function FerriteToComodo(grid)
    CellType = eltype(grid.cells)
    
    # Extract vertices (dimension depends on cell type)
    dim = CellType <: Union{Ferrite.Quadrilateral, Ferrite.Triangle} ? 2 : 3
    PointType = Point{dim, Float64}
    V = [PointType(Ferrite.get_node_coordinate(n)...) for n in grid.nodes]
    
    # Handle volumetric elements (3D)
    if CellType <: Union{Ferrite.Hexahedron, Ferrite.Tetrahedron}
        # Extract elements
        ElementType = CellType <: Ferrite.Hexahedron ? Hex8{Int64} : Tet4{Int64}
        E = [ElementType(cell.nodes) for cell in grid.cells]
        
        # Compute faces and boundary
        numElements = Ferrite.getncells(grid)
        CF_type = repeat(1:6, numElements)
        
        F = element2faces(E)
        F_uni, indUni, c_uni = gunique(F; 
                                        return_index=Val(true), 
                                        return_counts=Val(true), 
                                        sort_entries=true)
        
        Lb = isone.(c_uni)
        Fb = F_uni[Lb]
        CFb_type = CF_type[indUni][Lb]
        
        return E, V, F, Fb, CFb_type
    
    # Handle surface elements (2D)
    elseif CellType <: Ferrite.Quadrilateral
        F = [QuadFace{Int64}(cell.nodes) for cell in grid.cells]
        return F, V
    
    elseif CellType <: Ferrite.Triangle
        F = [TriangleFace{Int64}(cell.nodes) for cell in grid.cells]
        return F, V
    
    else
        error("Unsupported cell type: $(CellType). Expected Hexahedron, Tetrahedron, Quadrilateral, or Triangle.")
    end
end
###########################################################
###########################################################
function faceset_to_cellid_faceid(E::Vector{<: AbstractElement{N, T}}, faceIndices) where N where T 
    element_type = eltype(E)
    if element_type <: Tet4{T}
        nf = 4        
    elseif element_type <: Tet10{T}
        nf = 4                
    elseif element_type <: Tet15{T}
        nf = 4                
    elseif element_type <: Hex8{T}
        nf = 6        
    elseif element_type <: Penta6{T}
        nf = 5
    end
    face_cell_id = Vector{Int}(undef, length(faceIndices))
    face_id = Vector{Int}(undef, length(faceIndices))
    @inbounds for (indFace, i) in enumerate(faceIndices)
        face_cell_id[indFace] = ceil(Int,i/nf)
        face_id[indFace] = mod1(i,nf)
    end
    return face_cell_id, face_id
end