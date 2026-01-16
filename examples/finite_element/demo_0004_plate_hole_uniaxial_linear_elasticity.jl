using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite

GLMakie.closeall()

w = 12.0
h = 6.0
pointSpacing = 0.5
V1 = rectanglepoints(w, h, pointSpacing; dir=:acw)
r = 1.5
n = ceil(Int, (2 * pi * r) / pointSpacing)
V2 = circlepoints(r, n; dir=:acw)
VT = (V1, V2,)
R = ([1, 2],)
P = (pointSpacing)

F,V,C = regiontrimesh(VT,R,P)
grid = ComodoToFerrite(F, V, Ferrite.Triangle)
addfacetset!(grid, "left", x -> x[1] ≈ -6.0)
addfacetset!(grid, "right", x -> x[1] ≈ 6.0);


## plot the mesh
M = GeometryBasics.Mesh(V, F)
fig_mesh = Figure(size=(800, 800))
ax = Axis(fig_mesh[1, 1], aspect=DataAspect(), xlabel="X", ylabel="Y",title="Mesh with Boundary Conditions")
hp1 = poly!(ax, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=1, shading =  true, transparency = false)
display(GLMakie.Screen(), fig_mesh)

function create_values()
    dim, order = 2, 1
    ip = Ferrite.Lagrange{Ferrite.RefTriangle,order}()^dim
    qr = Ferrite.QuadratureRule{Ferrite.RefTriangle}(2)
    qr_face = Ferrite.FacetQuadratureRule{Ferrite.RefTriangle}(1)
    cell_values = Ferrite.CellValues(qr, ip)
    facet_values = Ferrite.FacetValues(qr_face, ip)
    return cell_values, facet_values
end
# Function to create DofHandler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefTriangle,1}()^2)
    Ferrite.close!(dh)
    return dh
end

function create_bc(dh, grid, displacement_prescribed)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0, 0.0], [1, 2]) # bcSupportList_XY
    add!(ch, dbc)
    dbc = Dirichlet(:u, getfacetset(grid, "right"), (x, t) -> [displacement_prescribed], [1]) # bcPrescribeList_X
    add!(ch, dbc)
    Ferrite.close!(ch)
    return ch
end

fig_bc = Figure(size=(800, 800))
ax2 = Axis(fig_bc[1, 1], aspect=DataAspect(), title="Boundary condition")
hp2 = poly!(ax2, M, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=2.0, shading=true, transparency=true)


facesset_bcSupportList_XY= get_boundary_points(grid, getfacetset(grid, "left"), Faces, Ferrite.Quadrilateral)
scatter!(ax2, facesset_bcSupportList_XY, color=:blue, markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label="bcSupportList_XY")

facesset_bcPrescribeList_X= get_boundary_points(grid, getfacetset(grid, "right"), Faces, Ferrite.Quadrilateral)
scatter!(ax2, facesset_bcPrescribeList_X, color=:green, markersize=15.0, marker=:circle, strokecolor=:black, strokewidth=2, label="bcPrescribeList_X")

axislegend(ax2, position=:rt, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(GLMakie.Screen(), fig_bc)

function get_material_matrix(E, ν)
    C_voigt = E * [1.0 ν 0.0; ν 1.0 0.0; 0.0 0.0 (1-ν)/2] / (1 - ν^2)
    return fromvoigt(SymmetricTensor{4,2}, C_voigt)
end

# function to assemble the local stiffness matrix for 2D Plane stress
function assemble_cell!(ke, cell_values, E, ν)
    C = get_material_matrix(E, ν)
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

# function to assemble the global stiffness matrix for 2D
function assemble_global!(K, dh, cell_values, E, ν)
    n_basefuncs = getnbasefunctions(cell_values)
    ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K)

    for cell in CellIterator(dh)
        reinit!(cell_values, cell)
        fill!(ke, 0.0)
        assemble_cell!(ke, cell_values, E, ν)
        assemble!(assembler, celldofs(cell), ke)
    end
    return K
end

function solveLinearElasticSteps(E, ν, grid, displacement_prescribed, numSteps)
    dh = create_dofhandler(grid)
    numNodes = length(grid.nodes)
    
    # Step 0: undeformed configuration        
    U0 = zeros(Point{2, Float64}, numNodes) # Initialise displacement vectors for initial state 
    U0_mag = zeros(Float64, numNodes) # Initial displacement magnitude for initial state 
    UT = [U0 for _ in 1:numSteps]   # Initialise displacement vectors for each step
    UT_mag = [U0_mag for _ in 1:numSteps] # Initialise displacement magnitudes for each step
    ut_mag_max = zeros(Float64,numSteps) # Max magnitude per step
    for i = 1:numSteps-1
        iStep = i
        println("Solving step $iStep of $numSteps")
        ch = create_bc(dh, grid, (iStep/(numSteps-1))*displacement_prescribed)

        
        cell_values, facet_values = create_values()

        K = allocate_matrix(dh)
        assemble_global!(K, dh, cell_values, E, ν);

        f_ext = zeros(ndofs(dh))
        apply!(K, f_ext, ch)

        # Solve linear system
        U = K \ f_ext

        # Current deformed configuration
        u_nodes = vec(evaluate_at_grid_nodes(dh, U, :u))
        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)

        disp_points = [Point{2, Float64}([ux[j], uy[j]]) for j in eachindex(ux)]
        UT[i+1] = disp_points
        UT_mag[i+1] = norm.(disp_points)
        ut_mag_max[i+1] = maximum(UT_mag[i+1])
    end
    return UT, UT_mag, ut_mag_max
end

sampleSize = 3.17
strainApplied = 0.05 # 
loadingOption = "tension" #  "tension" or "compression"

E = 500.0e3 # MPa
ν = 0.3 # Poisson's ratio

if loadingOption == "tension"
    displacement_prescribed = strainApplied * sampleSize
elseif loadingOption == "compression"
    displacement_prescribed = -strainApplied * sampleSize
end

numSteps = 20

UT, UT_mag, ut_mag_max = solveLinearElasticSteps(E, ν, grid, displacement_prescribed, numSteps)

# Create displaced mesh per step
scale = 5.0

VV = [Point{2,Float64}(e[1], e[2]) for e in V]
VT = [VV .+ scale .* UT[i] for i in 1:numSteps]

min_p = minp([minp(V) for V in VT])
max_p = maxp([maxp(V) for V in VT])

# === Visualization setup ===
fig_disp = Figure(size=(1000,600))
stepStart = 2  # Start at undeformed
ax3 = Axis(fig_disp[1, 1], aspect=DataAspect(), title = "Step: $stepStart")

xlims!(ax3, min_p[1], max_p[1])
ylims!(ax3, min_p[2], max_p[2])
hp = poly!(ax3, GeometryBasics.Mesh(VT[stepStart], F), 
               strokewidth = 2,
               color = UT_mag[stepStart], 
               transparency = false, 
               colormap = Reverse(:Spectral), 
               colorrange = (0, maximum(ut_mag_max)))

Colorbar(fig_disp[1, 2], hp.plots[1], label = "Displacement magnitude [mm]") 

incRange = 1:numSteps
hSlider = Slider(fig_disp[2, 1], range = incRange, startvalue = stepStart - 1, linewidth = 30)

on(hSlider.value) do stepIndex     
    hp[1] = GeometryBasics.Mesh(VT[stepIndex], F)
    hp.color = UT_mag[stepIndex]
    ax3.title = "Step: $stepIndex"
end

slidercontrol(hSlider, ax3)
display(GLMakie.Screen(), fig_disp)