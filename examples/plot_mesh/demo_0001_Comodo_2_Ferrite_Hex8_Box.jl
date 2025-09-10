using Revise
using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Ferrite
using ComodoFerrite
using Colors


pointSpacing = 2.0
boxDim = [2,2,2]
boxEl = [6,6,6] 
E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl) 
grid = ComodoToFerrite(E, V,Fb, CFb_type, Ferrite.Hexahedron)
addnodeset!(grid, "traction", x -> x[1] ≈ 1.0) # we add a nodeset to check the boundary condition

# Function to create cell and facet values
function create_values()
    order = 1
    dim = 3
    ip = Lagrange{RefHexahedron,order}()^dim
    qr = QuadratureRule{RefHexahedron}(2)
    qr_face = FacetQuadratureRule{RefHexahedron}(1)
    cell_values = CellValues(qr, ip)
    facet_values = FacetValues(qr_face, ip)
    return cell_values, facet_values
end

# Function to create a DOF handler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron,1}()^3)
    Ferrite.close!(dh)
    return dh
end

# Function to create boundary conditions
function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0, 0.0, 0.0], [1,2,3]) ## change it to other (right, top, etc)
    add!(ch, dbc)
    Ferrite.close!(ch)
    return ch
end
dh = create_dofhandler(grid)
ch = create_bc(dh, grid)

## Visualize mesh 
GLMakie.closeall()

fig = Figure(size = (1200,800))

ax1 = AxisGeom(fig[1, 1], title = "Hex8 mesh", azimuth = -0.1π, elevation = 0.1π)
hp1 = meshplot!(ax1, F, V,color=:gray,  strokecolor=:black, strokewidth=3.0, shading =  false, transparency = false)

ax2 = AxisGeom(fig[1, 2], title = "Boundary condition", azimuth = -0.1π, elevation = 0.1π)
hp2 = meshplot!(ax2, Fb, V, color=(Gray(0.95), 0.3),  strokecolor=:black, strokewidth=2.0, shading=true, transparency = true)
 
## Visualize boundary condition
### First convert the face_index and node_index to points
facesset = get_boundary_points(grid, getfacetset(grid, "left"), Faces, Ferrite.Hexahedron)
nodesset = get_boundary_points(grid, getnodeset(grid, "traction"), Nodes, Ferrite.Hexahedron)

scatter!(ax2, facesset, color=:blue,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Fixed XYZ")
scatter!(ax2, nodesset, color=:red,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Traction")
axislegend(ax2, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)

display(fig)

