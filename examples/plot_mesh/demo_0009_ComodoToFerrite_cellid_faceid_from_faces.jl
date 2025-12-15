using ComodoFerrite.Comodo
using ComodoFerrite.Comodo.GLMakie
using ComodoFerrite.Comodo.GLMakie.Colors
using ComodoFerrite.Comodo.GeometryBasics
using ComodoFerrite.Comodo.Statistics

using ComodoFerrite
using ComodoFerrite.Ferrite

GLMakie.closeall()

for testCase = 1:2
    if testCase == 1
        pointSpacing = 0.5
        boxDim = [2.5,3.1,4] # Dimensionsions for the box in each direction
        boxEl = ceil.(Int,boxDim./pointSpacing) # Number of elements to use in each direction 
        E, V, F, _, _ = hexbox(boxDim,boxEl)
        nf = 6 
    elseif testCase == 2
        boxDim = [2.5,3.1,4] # Dimensions for the box in each direction
        pointSpacing = 0.5
        E, V, _, _ = tetbox(boxDim,pointSpacing)
        F = element2faces(E)
        nf = 4 
    end

    # Get boundary face indices 
    indBoundaryFaces = boundaryfaceindices(F)

    Fb = F[indBoundaryFaces]

    ZF = [mean(V[f])[3] for f in Fb]
    zMax = maximum(ZF)
    searchTol = 1e-6
    B = [z>(zMax-searchTol) for z in ZF]
    indBoundaryFaces_top = indBoundaryFaces[B]
    Fb_top = Fb[B]

    face_cell_id, face_id = faceset_to_cellid_faceid(E, indBoundaryFaces_top)

    # Visualisation
    Fbs_top,Vbs_top = separate_vertices(Fb_top,V)
    Cbs_top_face_cell_id = simplex2vertexdata(Fbs_top, face_cell_id)
    Cbs_top_face_id = simplex2vertexdata(Fbs_top, face_id)

    fig = Figure(size=(1600,800))

    ax1 = AxisGeom(fig[1, 1], title = "Faces featuring element ids")
    hp1 = meshplot!(ax1, Fbs, Vbs; strokewidth=0.0, color=(:white, 0.25), transparency=true)
    hp2 = meshplot!(ax1, Fbs_top, Vbs_top; strokewidth=0.0, color=Cbs_top_face_cell_id, colorrange=(1,length(E)))
    Colorbar(fig[1, 2], hp2)    

    ax2 = AxisGeom(fig[1, 3], title = "Faces featuring face ids")
    hp3 = meshplot!(ax2, Fbs, Vbs; strokewidth=0.0, color=(:white, 0.25), transparency=true)
    hp4 = meshplot!(ax2, Fbs_top, Vbs_top; strokewidth=0.0, color=Cbs_top_face_id, colorrange=(1,nf))    
    Colorbar(fig[1, 4], hp4)

    screen = display(GLMakie.Screen(), fig)
    GLMakie.set_title!(screen, "testCase = $testCase")
end