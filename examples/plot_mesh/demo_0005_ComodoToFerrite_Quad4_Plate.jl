using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite

plateDim = [1.0, 1.0]
plateElem = [10, 10]
orientation = :up
F, V, Eb, Cb = quadplate(plateDim, plateElem; orientation = orientation)

grid = ComodoToFerrite(F, V)

Fb_left = Eb[Cb.==4]

addface!(grid, "left", Fb_left)
GLMakie.closeall()

M = GeometryBasics.Mesh(V, F)
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y",title="Mesh with Boundary Conditions")
hp1 = meshplot!(ax, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=1, shading =  true, transparency = false)

facesset = get_boundary_points(grid, getfacetset(grid, "left"))
scatter!(ax, facesset, color=:blue, markersize=15.0, marker=:circle, label = "Fixed XY")

display(GLMakie.Screen(), fig)


