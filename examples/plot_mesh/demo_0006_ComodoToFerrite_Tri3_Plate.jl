using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite

plateDim1 = [20.0,24.0]
pointSpacing1 = 2.0

orientation1 = :up
F, V, Eb, Cb = triplate(plateDim1, pointSpacing1; orientation= orientation1, return_boundary_edges =  Val(true) )


grid = ComodoToFerrite(F, V)


GLMakie.closeall()

M = GeometryBasics.Mesh(V, F)
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y",title="Mesh with Boundary Conditions")
poly!(ax, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=1, shading =  true, transparency = false)

facesset = get_boundary_points(grid, getfacetset(grid, "right"), Faces, Ferrite.Quadrilateral)
scatter!(ax, facesset, color=:blue, markersize=15.0, marker=:circle, label = "Fixed XY")

display(GLMakie.Screen(), fig)