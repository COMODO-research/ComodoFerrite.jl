using ComodoFerrite.Comodo
using ComodoFerrite.Comodo.GLMakie
using ComodoFerrite.Comodo.GLMakie.Colors
using ComodoFerrite.Comodo.GeometryBasics
using ComodoFerrite.Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite

boxDim = [2,2,2]
boxEl = [6,6,6] 
E, V, F, Fb, Cb = hexbox(boxDim,boxEl) 
grid = ComodoToFerrite(E, V, Ferrite.Hexahedron; Fb, Cb)
grid.facetsets
addnodeset!(grid, "traction", x -> x[1] ≈ 1.0) # we add a nodeset to check the boundary condition


## Visualize mesh 
GLMakie.closeall()

M = GeometryBasics.Mesh(V, F, normal = face_normals(V, F))
Mb = GeometryBasics.Mesh(V, Fb, normal = face_normals(V, Fb))

fig = Figure(size = (1200,800))

ax1 = AxisGeom(fig[1, 1], title = "Hex8 mesh", azimuth = -0.2π, elevation = 0.1π)
hp1 = meshplot!(ax1, M,color=:gray,  strokecolor=:black, strokewidth=3.0, shading =  true, transparency = false)

ax2 = AxisGeom(fig[1, 2], title = "Boundary condition", azimuth = -0.2π, elevation = 0.1π)
hp2 = meshplot!(ax2, Mb, color=(Gray(0.95), 0.3),  strokecolor=:black, strokewidth=2.0, shading=true, transparency = true)
 
## Visualize boundary condition
### First convert the face_index and node_index to points
facesset = get_boundary_points(grid, getfacetset(grid, "left"), Faces, Ferrite.Hexahedron)
nodesset = get_boundary_points(grid, getnodeset(grid, "traction"), Nodes, Ferrite.Hexahedron)

scatter!(ax2, facesset, color=:blue,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Fixed XYZ")
scatter!(ax2, nodesset, color=:red,markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "Traction")
axislegend(ax2, position=:rb, backgroundcolor=(:white, 0.7), framecolor=:gray)

display(fig)

