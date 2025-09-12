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


