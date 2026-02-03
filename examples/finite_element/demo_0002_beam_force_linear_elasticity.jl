using Comodo
using Comodo.GLMakie
using Comodo.GLMakie.Colors
using Comodo.GeometryBasics
using Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite
#=
Finite element for beam in linear elastic regime 
=#

## GLMakie setting 
GLMakie.closeall()
# GLMakie.activate!()
GLMakie.activate!()
## Mesh 
pointSpacing = 3.0
boxDim = [10.0, 40.0, 10.0]
boxEl = ceil.(Int64, boxDim ./ pointSpacing)
E, V, F, Fb, Cb = hexbox(boxDim, boxEl)
grid = ComodoToFerrite(E, V)

Fb_back = Fb[Cb.==4]   
addface!(grid , "back", Fb_back) 

Fb_front = Fb[Cb.==3]  
addface!(grid , "front", Fb_front) 

## plot the mesh
fig_mesh = Figure(size=(800, 800))
ax1 = AxisGeom(fig_mesh[1, 1], title="Hex8 mesh")
hp1 = meshplot!(ax1, F, V, color=:gray, strokecolor=:black, strokewidth=3.0, shading=false, transparency=false)
# xlims!(ax1, -8, 8)
# ylims!(ax1, -22, 22)
# zlims!(ax1, -8, 8)

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
    dbc = Dirichlet(:u, getfacetset(grid, "back"), (x, t) -> [0.0, 0.0, 0.0], [1, 2, 3]) # bcSupport
    add!(ch, dbc)
    Ferrite.close!(ch)
    return ch
end
## plot boundary condition
fig_bc = Figure(size=(800, 800))
ax2 = AxisGeom(fig_bc[1, 1], title="Boundary condition")
hp2 = meshplot!(ax2, Fb, V, color=(Gray(0.95), 0.3), strokecolor=:black, strokewidth=2.0, shading=true, transparency=true)
# xlims!(ax2, -8, 8)
# ylims!(ax2, -22, 22)
# zlims!(ax2, -8, 8)

facesset_bcSupport = get_boundary_points(grid, getfacetset(grid, "back"))
scatter!(ax2, facesset_bcSupport, color=:blue, markersize=15.0, marker=:circle,label="bcSupport")

facesset_bcPrescribe = get_boundary_points(grid, getfacetset(grid, "front"))
scatter!(ax2, facesset_bcPrescribe, color=:green, markersize=15.0, marker=:circle,label="bcPrescribe")

axislegend(ax2, position=:rt, backgroundcolor=(:white, 0.7), framecolor=:gray)
display(GLMakie.Screen(), fig_bc)

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

## Traction force
function assemble_external_forces!(f_ext, dh, facetset, facetvalues, prescribed_traction)
    # Create a temporary array for the facet's local contributions to the external force vector
    fe_ext = zeros(getnbasefunctions(facetvalues))
    for facet in FacetIterator(dh, facetset)
        # Update the facetvalues to the correct facet number
        reinit!(facetvalues, facet)
        # Reset the temporary array for the next facet
        fill!(fe_ext, 0.0)
        # Access the cell's coordinates
        for qp in 1:getnquadpoints(facetvalues)
            # Calculate the global coordinate of the quadrature point.
            # Get the integration weight for the current quadrature point.
            dΓ = getdetJdV(facetvalues, qp)
            for i in 1:getnbasefunctions(facetvalues)
                Nᵢ = shape_value(facetvalues, qp, i)
                fe_ext[i] += prescribed_traction ⋅ Nᵢ * dΓ
            end
        end
        # Add the local contributions to the correct indices in the global external force vector
        assemble!(f_ext, celldofs(facet), fe_ext)
    end
    return f_ext
end

function solveLinearElasticSteps(E, ν, grid, prescribed_traction, numSteps)
    dh = create_dofhandler(grid)
    numNodes = length(grid.nodes)
    
    # Step 0: undeformed configuration        
    U0 = zeros(Point{3, Float64}, numNodes) # Initialise displacement vectors for initial state 
    U0_mag = zeros(Float64, numNodes) # Initial displacement magnitude for initial state 
    UT = [U0 for _ in 1:numSteps]   # Initialise displacement vectors for each step
    UT_mag = [U0_mag for _ in 1:numSteps] # Initialise displacement magnitudes for each step
    ut_mag_max = zeros(Float64,numSteps) # Max magnitude per step
    for i = 1:numSteps-1
        iStep = i
        println("Solving step $iStep of $numSteps")
        ch = create_bc(dh, grid)

        cell_values, facet_values = create_values()

        K = allocate_matrix(dh)
        assemble_global_3d!(K, dh, cell_values, E, ν);

        f_ext = zeros(ndofs(dh))
        facetset = getfacetset(grid, "front")
        assemble_external_forces!(f_ext, dh, facetset, facet_values,  (iStep/(numSteps-1)) .*prescribed_traction)
        apply!(K, f_ext, ch)

        # Solve linear system
        U = K \ f_ext

        # Current deformed configuration
        u_nodes = vec(evaluate_at_grid_nodes(dh, U, :u))
        ux = getindex.(u_nodes, 1)
        uy = getindex.(u_nodes, 2)
        uz = getindex.(u_nodes, 3)

        disp_points = [Point{3, Float64}([ux[j], uy[j], uz[j]]) for j in eachindex(ux)]
        UT[i+1] = disp_points
        UT_mag[i+1] = norm.(disp_points)
        ut_mag_max[i+1] = maximum(UT_mag[i+1])
    end
    return UT, UT_mag, ut_mag_max
end

E = 500.0e3 # MPa
ν = 0.3 # Poisson's ratio

prescribed_traction = (0.0, 0.0, -30)
numSteps = 20

UT, UT_mag, ut_mag_max = solveLinearElasticSteps(E, ν, grid, prescribed_traction, numSteps)


# Create displaced mesh per step
scale = 10.0
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