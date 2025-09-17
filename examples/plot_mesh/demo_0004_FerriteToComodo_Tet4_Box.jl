using Revise , Comodo , Comodo.GLMakie ,Comodo.GeometryBasics , Comodo.Statistics
using Ferrite, ComodoFerrite , Colors

Nx = 5;  Ny = 5 ; Nz = 5
Lx = 2.0; Ly = 2.0; Lz = 2.0
left = Ferrite.Vec(0.0, 0.0, 0.0)
right = Ferrite.Vec(Lx, Ly, Lz)
grid = generate_grid(Ferrite.Tetrahedron, (Nx, Ny, Nz), left, right)

E , V, F, Fb, CFb_type   = FerriteToComodo(grid, Ferrite.Tetrahedron)
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