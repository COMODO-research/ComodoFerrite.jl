using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite

boxDim = [4.0 , 4.0 , 4.0] # Dimensions for the box in each direction
pointSpacing = 0.9

E, V, Fb, Cb = tetbox(boxDim, pointSpacing)

F = element2faces(E)

grid = ComodoToFerrite( E, V,  Ferrite.Tetrahedron ; Fb, Cb)
addfacetset!(grid, "traction", x -> x[3] ≈ 2.0)


## Visualize mesh 
GLMakie.closeall()

M = GeometryBasics.Mesh(V, F, normal = face_normals(V, F))
Mb = GeometryBasics.Mesh(V, Fb, normal = face_normals(V, Fb))

fig = Figure(size = (1200,800))

ax1 = AxisGeom(fig[1, 1], title = "Tet4 mesh", azimuth = -0.2π, elevation = 0.1π)
hp1 = meshplot!(ax1, M ,color=:gray,  strokecolor=:black, strokewidth=3.0, shading =  false, transparency = false)

ax2 = AxisGeom(fig[1, 2], title = "Boundary condition", azimuth = -0.2π, elevation = 0.1π)
hp2 = meshplot!(ax2, Mb, color=(Gray(0.95), 0.3),  strokecolor=:black, strokewidth=2.0, shading=true, transparency = true)
 
## Visualize boundary condition
### First convert the face_index and node_index to points
facesset = get_boundary_points(grid, getfacetset(grid, "traction"), Faces, Ferrite.Tetrahedron)
nodesset = get_boundary_points(grid, getfacetset(grid, "left"), Faces, Ferrite.Tetrahedron)

scatter!(ax2, facesset, color=:blue,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Traction")
scatter!(ax2, nodesset, color=:red,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Fixed_XYZ")
axislegend(ax2, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)

display(fig)


