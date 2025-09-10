using Revise , Comodo , Comodo.GLMakie ,Comodo.GeometryBasics , Comodo.Statistics
using Ferrite, ComodoFerrite , Colors

Nx = 5;  Ny = 5 ; Nz = 5
Lx = 2.0; Ly = 2.0; Lz = 2.0
left = Ferrite.Vec(0.0, 0.0, 0.0)
right = Ferrite.Vec(Lx, Ly, Lz)
grid = generate_grid(Ferrite. Hexahedron, (Nx, Ny, Nz), left, right)

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
    dbc = Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> [0.0, 0.0, 0.0], [1,2,3])
    dbc = Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> [0.0], [1])
    add!(ch, dbc)
    Ferrite.close!(ch)
    return ch
end
dh = create_dofhandler(grid)
ch = create_bc(dh, grid)


E , V, F, Fb  = FerriteToComodo(grid, Ferrite.Hexahedron)

## Visualize mesh 
GLMakie.closeall()

fig = Figure(size = (1200,800))

ax1 = AxisGeom(fig[1, 1], title = "Hex8 mesh", azimuth = -0.2π, elevation = 0.1π)
hp1 = meshplot!(ax1, F, V,color=:gray,  strokecolor=:black, strokewidth=3.0, shading =  false, transparency = false)

ax2 = AxisGeom(fig[1, 2], title = "Boundary condition", azimuth = -0.2π, elevation = 0.1π)
hp2 = meshplot!(ax2, Fb, V, color=(Gray(0.95), 0.3),  strokecolor=:black, strokewidth=2.0, shading=true, transparency = true)
 
facesset_bottom = get_boundary_points(grid, getfacetset(grid, "bottom"), Faces, Ferrite.Hexahedron)
facesset_top = get_boundary_points(grid, getfacetset(grid, "top"), Faces, Ferrite.Hexahedron)


scatter!(ax2, facesset_bottom, color=:blue,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Fixed XYZ")
scatter!(ax2, facesset_top, color=:red,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Fixed_x")
axislegend(ax2, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)

display(fig)