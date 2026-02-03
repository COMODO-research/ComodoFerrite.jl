using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite

Nx = 5;  Ny = 5 ; Nz = 5
Lx = 2.0; Ly = 2.0; Lz = 2.0
left = Ferrite.Vec(0.0, 0.0, 0.0)
right = Ferrite.Vec(Lx, Ly, Lz)
grid = generate_grid(Ferrite.Hexahedron, (Nx, Ny, Nz), left, right)

E , V, F, Fb, CFb_type   = FerriteToComodo(grid)

## Visualize mesh 
GLMakie.closeall()

M = GeometryBasics.Mesh(V, F, normal = face_normals(V, F))
Mb = GeometryBasics.Mesh(V, Fb, normal = face_normals(V, Fb))

fig = Figure(size = (1200,800))

ax1 = AxisGeom(fig[1, 1], title = "Hex8 mesh")
hp1 = meshplot!(ax1, M, color=:gray,  strokecolor=:black, strokewidth=3.0, shading =  false, transparency = false)

ax2 = AxisGeom(fig[1, 2], title = "Boundary condition")
hp2 = meshplot!(ax2, Mb, color=(Gray(0.95), 0.3),  strokecolor=:black, strokewidth=2.0, shading=true, transparency = true)
 
facesset_bottom = get_boundary_points(grid, getfacetset(grid, "bottom"))
facesset_top = get_boundary_points(grid, getfacetset(grid, "top"))


scatter!(ax2, facesset_bottom, color=:blue,markersize=15.0, marker=:circle,  label = "Fixed XYZ")
scatter!(ax2, facesset_top, color=:red,markersize=15.0, marker=:circle, label = "Fixed_x")
axislegend(ax2, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)

display(GLMakie.Screen(), fig)