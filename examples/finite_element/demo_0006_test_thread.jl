### check the Ferrite for Thread
using ComodoFerrite
using ComodoFerrite.Ferrite
using ComodoFerrite.Comodo
using ComodoFerrite.Comodo.GeometryBasics

using SparseArrays
using OhMyThreads
using Base.Threads

### first in the terminal:  julia -t auto or  julia -t NumberofThread
### before run check the number of thread by Threads.nthreads()
### then run the code in REPL by include("demo_0006_test_thread.jl")
### 
    

#### we can check the number of CPUs by Sys.CPU_THREADS


function create_grid()
    nx, ny  = (1000 , 1000)
    Lx, Ly  = (2. , 2.)
    corners = [
        ComodoFerrite.Ferrite.Vec{2}((0.0, 0.0)), ComodoFerrite.Ferrite.Vec{2}((Lx, 0.0)),
        ComodoFerrite.Ferrite.Vec{2}((Lx, Ly)), ComodoFerrite.Ferrite.Vec{2}((0.0, Ly))
    ]
    grid = ComodoFerrite.Ferrite.generate_grid(ComodoFerrite.Ferrite.Quadrilateral, (nx, ny), corners)
    colors = create_coloring(grid)

    return grid, colors
end

function create_material_stiffness()
    E = 10e9
    ν = 0.3
    C_voigt = E * [1.0 ν 0.0; ν 1.0 0.0; 0.0 0.0 (1-ν)/2] / (1 - ν^2)
    return fromvoigt(SymmetricTensor{4,2}, C_voigt)
end

function create_values()
    order = 1
    dim = 2
    ip = Lagrange{RefQuadrilateral,order}()^dim
    qr = QuadratureRule{RefQuadrilateral}(2)
    qr_face = FacetQuadratureRule{RefQuadrilateral}(1)
    cellvalues = CellValues(qr, ip)
    facetvalues = FacetValues(qr_face, ip)
    return cellvalues, facetvalues
end
function create_dofhandler(grid)
    dh = ComodoFerrite.Ferrite.DofHandler(grid)
    ComodoFerrite.Ferrite.add!(dh, :u, ComodoFerrite.Ferrite.Lagrange{ComodoFerrite.Ferrite.RefQuadrilateral,1}()^2)
    ComodoFerrite.Ferrite.close!(dh)
    return dh
end

function assemble_cell!(dh, Ke::Matrix, fe::Vector, cellvalues::CellValues, facetvalues::FacetValues,
    C::SymmetricTensor, facetset, traction::ComodoFerrite.Ferrite.Vec)
    fill!(Ke, 0)
    fill!(fe, 0)
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        for i in 1:getnbasefunctions(cellvalues)
            ∇δui = shape_symmetric_gradient(cellvalues, q_point, i)
            for j in 1:getnbasefunctions(cellvalues)
                ∇uj = shape_symmetric_gradient(cellvalues, q_point, j)
                Ke[i, j] += (∇δui ⊡ C ⊡ ∇uj) * dΩ
            end
        end
    end
    for facet in FacetIterator(dh, facetset)
        reinit!(facetvalues, facet)
        for qp in 1:getnquadpoints(facetvalues)
            dΓ = getdetJdV(facetvalues, qp)
            for i in 1:getnbasefunctions(facetvalues)
                Nᵢ = shape_value(facetvalues, qp, i)
                fe[i] += traction ⋅ Nᵢ * dΓ
            end
        end
    end
    return Ke, fe
end

struct ScratchData{CC, CV, FV, FS, T, A}
    cell_cache::CC
    cellvalues::CV
    facetvalues::FV
    facetset::FS
    Ke::Matrix{T}
    fe::Vector{T}
    assembler::A
end

function ScratchData(dh::DofHandler, K::SparseMatrixCSC, f::Vector,
                     cellvalues::CellValues, facetvalues::FacetValues, facetset)
    cell_cache = CellCache(dh)
    n = ndofs_per_cell(dh)
    T = eltype(f)
    Ke = zeros(T, n, n)
    fe = zeros(T, n)
    asm = start_assemble(K, f; fillzero = false)
    ScratchData(cell_cache, copy(cellvalues), copy(facetvalues), copy(facetset), Ke, fe, asm)
end

function assemble_global!(
        K::SparseMatrixCSC, f::Vector, dh::DofHandler, colors,
        cellvalues_template::CellValues, facetvalues_template, facetset_template; ntasks = Threads.nthreads()
    )
    # Zero-out existing data in K and f
    _ = start_assemble(K, f)
    # Body force and material stiffness
    C = create_material_stiffness()
    traction = ComodoFerrite.Ferrite.Vec(0.0, 1.0)
    # Loop over the colors
    for color in colors
        # Dynamic scheduler spawning `ntasks` tasks where each task will process a chunk of
        # (roughly) equal number of cells (`length(color) ÷ ntasks`).
        scheduler = OhMyThreads.DynamicScheduler(; ntasks)
        # Parallelize the loop over the cells in this color
        OhMyThreads.@tasks for cellidx in color
            # Tell the @tasks loop to use the scheduler defined above
            @set scheduler = scheduler
            # Obtain a task local scratch and unpack it
            @local scratch = ScratchData(dh, K, f, cellvalues_template, facetvalues_template, facetset_template )
            (; cell_cache, cellvalues, facetvalues, facetset, Ke, fe, assembler) = scratch
            # Reinitialize the cell cache and then the cellvalues
            reinit!(cell_cache, cellidx)
            reinit!(cellvalues, cell_cache)
            # Compute the local contribution of the cell
            assemble_cell!(dh , Ke, fe, cellvalues, facetvalues, C, facetset, traction)
            # Assemble local contribution
            assemble!(assembler, celldofs(cell_cache), Ke, fe)
        end
    end
    return K, f
end

function main(; ntasks = Threads.nthreads())
    # Interpolation, quadrature and cellvalues
    # Grid, colors and DofHandler
    grid, colors = create_grid()
    dh = create_dofhandler(grid)
    cellvalues, facetvalues = create_values()
    facetset = getfacetset(grid, "right")

    # Global matrix and vector
    K = allocate_matrix(dh)
    f = zeros(ndofs(dh))
    # Compile it
    assemble_global!( K, f, dh, colors, cellvalues, facetvalues, facetset; ntasks = ntasks)
    # Time it
    @time assemble_global!( K, f, dh, colors,cellvalues, facetvalues, facetset; ntasks = ntasks)
    return
end
##### test different ntasks for a thread number like 4
main(; ntasks = 1) 
main(; ntasks = 2) 
main(; ntasks = 3) 
main(; ntasks = 4) 