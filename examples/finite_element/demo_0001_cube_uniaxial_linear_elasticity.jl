using Revise
using Comodo
using Comodo.GLMakie
using Comodo.GeometryBasics
using Comodo.Statistics
using Ferrite
using ComodoFerrite
using Colors
#=
Finite element for cube uniaxial in linear elastic regime 
=#

## GLMakie setting 
GLMakie.closeall()
GLMakie.activate!( fullscreen=true)

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
function create_bc(dh, grid, displacement_prescribed)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(grid, "bottom"), (x, t) -> [0.0], [3]) # bcSupportList_Z
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "front"), (x, t) -> [0.0], [2]) # bcSupportList_Y
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0], [1]) # bcSupportList_X
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "top"), (x, t) -> [displacement_prescribed], [3]) # bcPrescribeList_Z
    add!(ch, dbc)
    Ferrite.close!(ch)
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

######################################################

## Finite element solver
# Define the 3D material stiffness matrix in Voigt notation
function get_material_matrix_3d(E, ν)
    # Define the 3D material stiffness matrix in Voigt notation
    C_voigt = E / ((1 + ν) * (1 - 2 * ν)) * [
        1-ν ν ν 0 0 0;
        ν 1-ν ν 0 0 0;
        ν ν 1-ν 0 0 0;
        0 0 0 (1-2*ν)/2 0 0;
        0 0 0 0 (1-2*ν)/2 0;
        0 0 0 0 0 (1-2*ν)/2
    ]
    return fromvoigt(SymmetricTensor{4,3}, C_voigt)
end

# function to assemble the local stiffness matrix for 3D
function assemble_cell_3d!(ke, cell_values, E, ν)
    C = get_material_matrix_3d(E, ν)
    for qp in 1:getnquadpoints(cell_values)
        dΩ = getdetJdV(cell_values, qp)
        for i in 1:getnbasefunctions(cell_values)
            ∇Ni = shape_gradient(cell_values, qp, i)
            for j in 1:getnbasefunctions(cell_values)
                ∇δNj = shape_symmetric_gradient(cell_values, qp, j)
                ke[i, j] += (∇Ni ⊡ C ⊡ ∇δNj) * dΩ
            end
        end
    end
    return ke
end
# function to assemble the global stiffness matrix for 3D
function assemble_global_3d!(K, dh, cell_values, E, ν)
    n_basefuncs = getnbasefunctions(cell_values)
    ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K)

    for cell in CellIterator(dh)
        reinit!(cell_values, cell)
        fill!(ke, 0.0)
        assemble_cell_3d!(ke, cell_values, E, ν)
        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end


sampleSize = 10.0
strainApplied = .5 # Equivalent linear strain
loadingOption = "tension" # "tension" or "compression"

if loadingOption == "tension"
    displacement_prescribed = strainApplied * sampleSize
elseif loadingOption == "compression"
    displacement_prescribed = -strainApplied * sampleSize
end

dh = create_dofhandler(grid)
ch = create_bc(dh, grid, displacement_prescribed)

cell_values, facet_values = create_values()
E = 500.0e3 # MPa
ν = 0.3 # Poisson's ratio

K = allocate_matrix(dh)
assemble_global_3d!(K, dh, cell_values, E, ν);


f_ext = zeros(ndofs(dh))
apply!(K, f_ext, ch)

# Solve linear system
U = K \ f_ext

# === Prepare data for single displacement step ===
totalSteps = 2  # Step 0 (undeformed) and Step 1 (deformed)

UT = Vector{Vector{Point{3, Float64}}}(undef, totalSteps)   # Displaced geometry
UT_mag = Vector{Vector{Float64}}(undef, totalSteps)         # Magnitude per node
ut_mag_max = zeros(totalSteps)                              # Max magnitude per step

# Step 0: undeformed configuration
UT[1] = [Point{3, Float64}([0.0, 0.0, 0.0]) for _ in V]
UT_mag[1] = zeros(length(V))
ut_mag_max[1] = 0.0

# Step 1: deformed configuration
u_nodes = vec(evaluate_at_grid_nodes(dh, U, :u))
ux = getindex.(u_nodes, 1)
uy = getindex.(u_nodes, 2)
uz = getindex.(u_nodes, 3)

disp_points = [Point{3, Float64}([ux[j], uy[j], uz[j]]) for j in eachindex(ux)]
UT[2] = disp_points

UT_mag[2] = norm.(disp_points)
ut_mag_max[2] = maximum(UT_mag[2])

# Create displaced mesh per step
scale = 1.0
VT = [V .+ scale .* UT[i] for i in 1:totalSteps]

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

incRange = 0:(totalSteps - 1)
hSlider = Slider(fig_disp[2, 1], range = incRange, startvalue = stepStart - 1, linewidth = 30)

on(hSlider.value) do stepIndex 
    step = stepIndex + 1  
    hp[1] = GeometryBasics.Mesh(VT[step], Fb)
    hp.color = UT_mag[step]
    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax3)
display(GLMakie.Screen(), fig_disp)

