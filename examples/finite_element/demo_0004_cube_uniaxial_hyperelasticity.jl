using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite


using FerriteHyperelastic  

const Faces = ComodoFerrite.Faces
const Nodes = ComodoFerrite.Nodes


# create the structure for saving fem input Similar to MATLAB
input = InputStruct()

## GLMakie setting 
GLMakie.closeall()


## Mesh 
boxDim = [10, 10, 10]
boxEl = [5, 5, 5]
E, V, F, Fb, Cb = hexbox(boxDim, boxEl)
grid = ComodoToFerrite(E, V, Ferrite.Hexahedron; Fb, Cb )

## plot the mesh
fig_mesh = Figure(size=(800, 800))
ax1 = AxisGeom(fig_mesh[1, 1], title="Hex8 mesh")
hp1 = meshplot!(ax1, F, V, color=:gray, strokecolor=:black, strokewidth=3.0, shading=false, transparency=false)
xlims!(ax1, -6, 6)
ylims!(ax1, -6, 6)
zlims!(ax1, -6, 6)

display(GLMakie.Screen(), fig_mesh)

## FEM Values (Interpolation and Quadrature points)
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


## Dof handler 
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron,1}()^3)
    Ferrite.close!(dh)
    return dh
end


## Dirichlet Boundary Condition
function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> [0.0], [3]) # bcSupportList_Z
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "front"), (x, t) -> [0.0], [2]) # bcSupportList_Y
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0], [1]) # bcSupportList_X
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> [t], [3]) # bcPrescribeList_Z
    add!(ch, dbc)
    Ferrite.close!(ch)
    Ferrite.update!(ch, 0.0)
    return ch
end
## plot boundary condition

fig_bc = Figure(size=(800, 800))
ax2 = AxisGeom(fig_bc[1, 1], title="Boundary condition")
hp2 = meshplot!(ax2, Fb, V, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=2.0, shading=true, transparency=true)
xlims!(ax2, -6, 6)
ylims!(ax2, -6, 6)
zlims!(ax2, -6, 6)

facesset_bcSupportList_Z = get_boundary_points(grid, getfacetset(grid, "bottom"), Faces, Ferrite.Hexahedron)
scatter!(ax2, facesset_bcSupportList_Z, color=:blue, markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label="bcSupportList_Z")

facesset_bcSupportList_Y = get_boundary_points(grid, getfacetset(grid, "front"), Faces, Ferrite.Hexahedron)
scatter!(ax2, facesset_bcSupportList_Y, color=:green, markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label="bcSupportList_Y")

facesset_bcSupportList_X = get_boundary_points(grid, getfacetset(grid, "left"), Faces, Ferrite.Hexahedron)
scatter!(ax2, facesset_bcSupportList_X, color=:red, markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label="bcSupportList_X")

facesset_bcPrescribeList_Z = get_boundary_points(grid, getfacetset(grid, "top"), Faces, Ferrite.Hexahedron)
scatter!(ax2, facesset_bcPrescribeList_Z, color=:black, markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label="bcPrescribeList_Z")
axislegend(ax2, position=:rt, backgroundcolor=(:white, 0.7), framecolor=:gray)


display(GLMakie.Screen(), fig_bc)

###### FEM 
function Ψ(C, μ, λ)
    J = sqrt(det(C))
    I1 = tr(C)
    
    return μ/2*(I1- 1) -μ*log(J) + λ/2*(log(J))^2
end

function constitutive_driver(C, μ, λ)
    ∂²Ψ∂C², ∂Ψ∂C = Tensors.hessian(y -> Ψ(y, μ, λ), C, :all)
    S = 2.0 * ∂Ψ∂C
    ∂S∂C = 2.0 * ∂²Ψ∂C²
    return S, ∂S∂C
end

function make_constitutive_driver(μ, λ)
    return C -> constitutive_driver(C, μ, λ)
end

input.model_type = :threeD   # or :plane_strain or :threeD
input.load_type = :displacement

input.E , input.ν = 1.0, 0.4
E = input.E
ν = input.ν
μ = E / (2 * (1 + ν))
λ = ν*E / ((1+ν)*(1-2ν))
input.material = make_constitutive_driver(μ, λ)
input.grid = grid
input.dh = create_dofhandler(grid)
input.ch = create_bc(input.dh, input.grid)
# Create CellValues and FacetValues
input.cell_values, input.facet_values = create_values()


input.dof_F = []
input.dof_U = []

sampleSize = 10.0
strainApplied = 0.5 # Equivalent linear strain
loadingOption = "compression" # "tension" or "compression"

if loadingOption == "tension"
    input.displacement = strainApplied * sampleSize
elseif loadingOption == "compression"
    input.displacement = -strainApplied * sampleSize
end

input.tol = 1e-6

input.filename = "2D_Hyper"
input.output_dir= "/Users/aminalibakhshi/Desktop/vtu_geo/"


maxIterPerInc,totalTime,initInc,minInc,maxInc,totalInc = initialize_solver(500,1.0,1e-2,1e-8,0.1,1000)
input.maxIterPerInc = maxIterPerInc
input.totalTime = totalTime
input.initInc = initInc
input.minInc = minInc
input.maxInc = maxInc
input.totalInc = totalInc

sol = run_fem(input)


incRange = length(sol.U_steps)
numSteps = incRange + 1  # Step 0 + displacement step
UT = Vector{Vector{Point{3, Float64}}}(undef, numSteps)         # Displaced geometry
UT_mag = Vector{Vector{Float64}}(undef, numSteps)               # Magnitude per node
ut_mag_max = zeros(numSteps)                                    # Max magnitude per step

# Step 0: undeformed configuration
UT[1] = [Point{3, Float64}([0.0, 0.0, 0.0]) for _ in V]  # zero displacement
UT_mag[1] = zeros(length(V))                       # zero magnitude
ut_mag_max[1] = 0.0

# Steps 1 to incRange: displacement steps
@inbounds for i in 1:incRange
    U = sol.U_steps[i]
    u_nodes = vec(evaluate_at_grid_nodes(input.dh, U, :u))
    ux = getindex.(u_nodes, 1)
    uy = getindex.(u_nodes, 2)
    uz = getindex.(u_nodes, 3)
    disp_points = [Point{3, Float64}([ux[j], uy[j], uz[j]]) for j in eachindex(ux)]
    UT[i+1] = disp_points
    UT_mag[i+1] = norm.(disp_points)
    ut_mag_max[i+1] = maximum(UT_mag[i+1])
end

scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:numSteps]

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

# === Visualization setup ===
fig_disp = Figure(size=(800, 800))
stepStart = 2  # Start at undeformed

ax3 = AxisGeom(fig_disp[1, 1], title = "Step: $stepStart", limits=(min_p[1], max_p[1], min_p[2], max_p[2], min_p[3], max_p[3]))
hp = meshplot!(ax3, Fb, VT[stepStart]; 
               strokewidth = 2,
               color = UT_mag[stepStart], 
               transparency = false, 
               colormap = Reverse(:Spectral), 
               colorrange = (0, maximum(ut_mag_max)))
Colorbar(fig_disp[1, 2], hp.plots[1], label = "Displacement magnitude [mm]") 

incRange = 1:numSteps
hSlider = Slider(fig_disp[2, 1], range = incRange, startvalue = stepStart - 1, linewidth = 30)

on(hSlider.value) do stepIndex     
    hp[1] = GeometryBasics.Mesh(VT[stepIndex], Fb)
    hp.color = UT_mag[stepIndex]
    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax3)
display(GLMakie.Screen(), fig_disp)

@info "max" minimum(minimum(UT[end]))
@info "max" maximum(maximum(UT[end]))