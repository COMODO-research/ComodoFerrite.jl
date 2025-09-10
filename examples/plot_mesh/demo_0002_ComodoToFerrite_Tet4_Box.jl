using Revise , Comodo , Comodo.GLMakie , Comodo.GeometryBasics , Comodo.Statistics
using Ferrite, ComodoFerrite , Colors
### The generated grid lacks the facetsets for the boundaries, so we add them by using Ferrite's addfacetset!.
boxDim = [2.5,3.1,4] # Dimensions for the box in each direction
pointSpacing = 0.5

E, V, Fb, Cb = tetbox(boxDim,pointSpacing)


Fbs, Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs,Cb)

grid = ComodoToFerrite(E, V,  Ferrite.Tetrahedron)
addnodeset!(grid, "Fixed_XYZ", x -> x[3] ≈ -2.0) # we add a nodeset to check the boundary condition
addfacetset!(grid, "traction", x -> x[3] ≈ 2.0)
# Function to create cell and facet values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefTetrahedron,order}()^dim
    qr = QuadratureRule{RefTetrahedron}(2)
    qr_face = FacetQuadratureRule{RefTetrahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

# Function to create a DOF handler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefTetrahedron,1}()^3)
    Ferrite.close!(dh)
    return dh
end

# Function to create boundary conditions
function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getnodeset(grid, "Fixed_XYZ"), (x, t) -> [0.0, 0.0, 0.0], [1,2,3]) ## change it to other (right, top, etc)
    add!(ch, dbc)
    Ferrite.close!(ch)
    return ch
end
dh = create_dofhandler(grid)
ch = create_bc(dh, grid)

## Visualize mesh 
GLMakie.closeall()

fig = Figure(size = (1200,800))

ax1 = AxisGeom(fig[1, 1], title = "Tet4 mesh", azimuth = -0.2π, elevation = 0.1π)
hp1 = meshplot!(ax1, Fbs, Vbs,color=:gray,  strokecolor=:black, strokewidth=3.0, shading =  false, transparency = false)

ax2 = AxisGeom(fig[1, 2], title = "Boundary condition", azimuth = -0.2π, elevation = 0.1π)
hp2 = meshplot!(ax2, Fb, V, color=(Gray(0.95), 0.3),  strokecolor=:black, strokewidth=2.0, shading=true, transparency = true)
 
## Visualize boundary condition
### First convert the face_index and node_index to points
facesset = get_boundary_points(grid, getfacetset(grid, "traction"), Faces, Ferrite.Tetrahedron)
nodesset = get_boundary_points(grid, getnodeset(grid, "Fixed_XYZ"), Nodes, Ferrite.Tetrahedron)

scatter!(ax2, facesset, color=:blue,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Traction")
scatter!(ax2, nodesset, color=:red,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Fixed_XYZ")
axislegend(ax2, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)

display(fig)


