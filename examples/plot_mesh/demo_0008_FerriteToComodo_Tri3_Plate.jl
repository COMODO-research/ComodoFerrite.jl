using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite

Nx, Ny  = (5 , 5)
Lx, Ly  = (2. , 2.)
left = Ferrite.Vec(0.0, 0.0)
right = Ferrite.Vec(Lx, Ly)
grid = generate_grid(Ferrite.Triangle, (Nx, Ny), left, right)
addnodeset!(grid, "traction", x -> x[1] ≈ 2.0) # we add a nodeset to check the boundary condition

F, V   = FerriteToComodo(grid, Ferrite.Triangle)

GLMakie.closeall()

M = GeometryBasics.Mesh(V, F)
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y",title="Mesh with Boundary Conditions")
poly!(ax, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=1, shading =  true, transparency = false)

facesset = get_boundary_points(grid, getnodeset(grid, "traction"), Nodes, Ferrite.Quadrilateral)
scatter!(ax, facesset, color=:blue, markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label = "traction")

display(GLMakie.Screen(), fig)