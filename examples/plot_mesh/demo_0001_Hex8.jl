using Revise
using Comodo
using Ferrite
using ComodoFerrite
pointSpacing = 2.0
boxDim = boxEl = [2,2,2]
boxEl = boxDim = [2,2,2] 
E,V,F,Fb,CFb_type = hexbox(boxDim,boxEl) 
grid = ComodoToFerrite(E, V,Fb, CFb_type,Ferrite.Hexahedron)

## test the faces using Ferite.jl functions
# Function to create cell and facet values
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

# Function to create a DOF handler
function create_dofhandler(grid)
    dh = Ferrite.DofHandler(grid)
    Ferrite.add!(dh, :u, Ferrite.Lagrange{Ferrite.RefHexahedron,1}()^3)
    Ferrite.close!(dh)
    return dh
end

# Function to create boundary conditions
function create_bc(dh, grid)
    ch = Ferrite.ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfacetset(grid, "left"), (x, t) -> [0.0, 0.0, 0.0], [1,2,3]) ## change it to other (right, top, etc)
    add!(ch, dbc)
    Ferrite.close!(ch)
    return ch
end
dh = create_dofhandler(grid)
ch = create_bc(dh, grid)

output_dir = "/Users/aminalibakhshi/Desktop/vtu_geo/"  # change to your direction
VTKGridFile(joinpath(output_dir, "boundary-conditions"), dh) do vtk
    Ferrite.write_constraints(vtk, ch)
end