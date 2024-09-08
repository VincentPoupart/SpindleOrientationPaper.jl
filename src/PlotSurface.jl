
function PlotSurfaceSpindle(gonad::String, node::Int, cell::String, wrlfiles)
    ####################### 
    #example of variable
    #=
    gonad = "r_2022-02-28_GSC_L4-06"
    node = 1 #please enter the node/timepoint here
    cell= " Cell_1"
    =#
    #######################
    fazer = string(gonad, cell)
    wrlbasenamelist = RemoveExt.(basename.(wrlfiles), ".")
    foo = gonad .== wrlbasenamelist
    wrl = unique(wrlfiles[foo])
    IndexedFaceSets = ExtractDataVRML(wrl[1])


    LengthPoints = length(IndexedFaceSets[node]["verts"])
    LengthCoordInd = length(IndexedFaceSets[node]["coordIndex"])


    X = Array{Float64}(undef, LengthPoints)
    Y = Array{Float64}(undef, LengthPoints)
    Z = Array{Float64}(undef, LengthPoints)
    I = Array{Int64}(undef, LengthCoordInd)
    J = Array{Int64}(undef, LengthCoordInd)
    K = Array{Int64}(undef, LengthCoordInd)

    for x in 1:length(IndexedFaceSets[node]["verts"])
        X[x] = IndexedFaceSets[node]["verts"][x][1]
        Y[x] = IndexedFaceSets[node]["verts"][x][2]
        Z[x] = IndexedFaceSets[node]["verts"][x][3]
    end
    for m in 1:length(IndexedFaceSets[node]["coordIndex"])
        I[m] = IndexedFaceSets[node]["coordIndex"][m][1] - 1#I,J,K must be 0 based python like
        J[m] = IndexedFaceSets[node]["coordIndex"][m][2] - 1
        K[m] = IndexedFaceSets[node]["coordIndex"][m][3] - 1
    end
    #X = convert(Array{Float64,1}, X)
    #Y = convert(Array{Float64,1}, Y)
    #Z = convert(Array{Float64,1}, Z)
    #I = convert(Array{Int64,1}, I)
    #J = convert(Array{Int64,1}, J)
    #K = convert(Array{Int64,1}, K)

    Xmax = maximum(X)
    Xmin = minimum(X)
    deltaX = Xmax - Xmin

    Ymax = maximum(Y)
    Ymin = minimum(Y)
    deltaY = Ymax - Ymin

    Zmax = maximum(Z)
    Zmin = minimum(Z)
    deltaZ = Zmax - Zmin

    X1 = filter(:FRAME => ==(node), Celloutput[fazer])
    spindle_mid_point = X1[1, :Spindle_Midpoint]
    cent = IndexedFaceSets[node]["centroids"]
    d = Vector{Float64}(undef, 0)
    for m in eachindex(cent)
        push!(d, Distance3D(cent[m], spindle_mid_point))
    end
    Ioo = sortperm(d)
    Sorted_Normals = IndexedFaceSets[node]["normals"][Ioo]
    Sorted_Areas = IndexedFaceSets[node]["areas"][Ioo]
    CumSumArea = cumsum(Sorted_Areas)
    distance = 1
    area = 0
    IndexDistance = 0
    if last(CumSumArea) > 100
        while area < 50
            global IndexDistance = findfirst(x -> x >= distance, d[Ioo])
            global area = CumSumArea[IndexDistance]
            global distance += 0.1
        end
        if IndexDistance !== nothing && d[Ioo][IndexDistance] <= 8
            spindle = [X1[1, :POSITION_X] - X1[1, :POSITION_X_1]; X1[1, :POSITION_Y] - X1[1, :POSITION_Y_1]; X1[1, :POSITION_Z] - X1[1, :POSITION_Z_1]]
            normal = sum(Sorted_Normals[1:IndexDistance])
            PatchCentroids = cent[Ioo][1:IndexDistance]
            PatchCoordIndex = IndexedFaceSets[node]["coordIndex"][Ioo][1:IndexDistance]
        end
    end

    X1A = Array{Float64}(undef, 2)
    Y1A = Array{Float64}(undef, 2)
    Z1A = Array{Float64}(undef, 2)

    Couleur = [1, 1]

    X1A[1] = X1[!, :POSITION_X][1]
    X1A[2] = X1[!, :POSITION_X_1][1]

    Y1A[1] = X1[!, :POSITION_Y][1]
    Y1A[2] = X1[!, :POSITION_Y_1][1]

    Z1A[1] = X1[!, :POSITION_Z][1]
    Z1A[2] = X1[!, :POSITION_Z_1][1]

    LengthCoordInd_patch = length(PatchCoordIndex)

    L = Array{Int64}(undef, LengthCoordInd_patch)
    M = Array{Int64}(undef, LengthCoordInd_patch)
    N = Array{Int64}(undef, LengthCoordInd_patch)

    for m in 1:LengthCoordInd_patch
        L[m] = PatchCoordIndex[m][1] .- 1
        M[m] = PatchCoordIndex[m][2] .- 1
        N[m] = PatchCoordIndex[m][3] .- 1
    end

    SMP = X1[!, :Spindle_Midpoint][1]
    SMP = SMP'
    AntiNormal = (-1 * normal) ./ 2

    E = Array{Float64}(undef, 2)
    F = Array{Float64}(undef, 2)
    G = Array{Float64}(undef, 2)

    E[1] = SMP[1]
    E[2] = SMP[1] + normal[1]
    F[1] = SMP[2]
    F[2] = SMP[2] + normal[2]
    G[1] = SMP[3]
    G[2] = SMP[3] + normal[3]

    EE = Array{Float64}(undef, 2)
    FF = Array{Float64}(undef, 2)
    GG = Array{Float64}(undef, 2)

    EE[1] = SMP[1]
    EE[2] = SMP[1] + AntiNormal[1]
    FF[1] = SMP[2]
    FF[2] = SMP[2] + AntiNormal[2]
    GG[1] = SMP[3]
    GG[2] = SMP[3] + AntiNormal[3]


    Xaxis = [1, 0, 0]
    Yaxis = [0, 1, 0]
    Zaxis = [0, 0, 1]
    MMM = copy(normal) #the normal to the plane
    xyz0 = SMP #% point inside the plane
    sz = 4 # size of the rectangle
    n = 1 #number of random points
    Q = nullspace(MMM')
    xyz = xyz0 + Q * ((rand(2, n) .- 0.5) * sz)
    SpinVect = ((xyz - xyz0) / (norm((xyz - xyz0))))'
    ortho1 = SpinVect[1, 1:3]
    ortho2 = cross(ortho1, normal) / norm(cross(ortho1, normal))
    a, b, c = 8, 0.5, 4  #ellipsoid principal axes
    # Euler angles
    α = 90 * pi / 180
    β = CalculerAngle(Yaxis, normal) - 90 * pi / 180
    γ = CalculerAngle(Xaxis, normal)
    ellipsoidM(u, v) = [a / 2 * cos(u) * sin(v), b / 2 * sin(u) * sin(v), c / 2 * cos(v)]    # ellipsoid
    RM(u, v) = RotXYZ(α, β, γ) * ellipsoidM(u, v) .+ C     # rotated ellipsoid
    a, b, c = 0.5, 4, 4 #ellipsoid principal axes
    lengthNormal = norm(normal)


    LenghtTranslation = 4
    f = LenghtTranslation / lengthNormal
    TranslationVect = normal * f
    Cent = spindle_mid_point[:]
    fazel = Cent + TranslationVect
    C = [fazel[1]; fazel[2]; fazel[3]] # center of ellipsoid
    u, v = range(0, 2π, length=72), range(0, π, length=72)
    xs, ys, zs = [[p[i] for p in RM.(u, v')] for i in 1:3]
    plot7 = PlotlyJS.surface(x=xs, y=ys, z=zs, opacity=0.1)

    LenghtTranslation = 2
    f = LenghtTranslation / lengthNormal
    TranslationVect = normal * f
    fazel = Cent + TranslationVect
    C = [fazel[1]; fazel[2]; fazel[3]] # center of ellipsoid
    u, v = range(0, 2π, length=72), range(0, π, length=72)
    xs, ys, zs = [[p[i] for p in RM.(u, v')] for i in 1:3]
    plot8 = PlotlyJS.surface(x=xs, y=ys, z=zs, opacity=0.2)

    LenghtTranslation = 0
    f = LenghtTranslation / lengthNormal
    TranslationVect = normal * f
    fazel = Cent + TranslationVect
    C = [fazel[1]; fazel[2]; fazel[3]] # center of ellipsoid
    u, v = range(0, 2π, length=72), range(0, π, length=72)
    xs, ys, zs = [[p[i] for p in RM.(u, v')] for i in 1:3]
    plot9 = PlotlyJS.surface(x=xs, y=ys, z=zs, opacity=1)

    LenghtTranslation = -2
    f = LenghtTranslation / lengthNormal
    TranslationVect = normal * f
    fazel = Cent + TranslationVect
    C = [fazel[1]; fazel[2]; fazel[3]] # center of ellipsoid
    u, v = range(0, 2π, length=72), range(0, π, length=72)
    xs, ys, zs = [[p[i] for p in RM.(u, v')] for i in 1:3]
    plot10 = PlotlyJS.surface(x=xs, y=ys, z=zs, opacity=0.5)

    LenghtTranslation = -4
    f = LenghtTranslation / lengthNormal
    TranslationVect = normal * f
    fazel = Cent + TranslationVect
    C = [fazel[1]; fazel[2]; fazel[3]] # center of ellipsoid
    u, v = range(0, 2π, length=72), range(0, π, length=72)
    xs, ys, zs = [[p[i] for p in RM.(u, v')] for i in 1:3]
    plot11 = PlotlyJS.surface(x=xs, y=ys, z=zs, opacity=0.5)

    camera = attr(eye=attr(x=1, y=5, z=1))

    layout3_5 = Layout(
        #width = 1000, height=800,
        title_text=gonad,
        title_x=0.5,
        titlefont_size="18",
        scene_aspectratio=attr(x=1, y=deltaY / deltaX, z=deltaZ / deltaX),
        #scene_aspectratio = :equal,
        #aspect_ratio = :equal,
        #scale = :none
    )



    plot1 = PlotlyJS.scatter(x=X1A, y=Y1A, z=Z1A, type="scatter3d", mode="markers", marker=attr(size=2, color=["red"; "red"], opacity=1), line=attr(color="red", width=2)) #Spindle poles
    plot2 = PlotlyJS.mesh3d(x=X, y=Y, z=Z, i=I, j=J, k=K, opacity=0.50, color="cyan")# mesh3d of the surface
    plot3 = PlotlyJS.mesh3d(x=X, y=Y, z=Z, i=L, j=M, k=N, opacity=0.5, color="blue")# plot patch
    plot4 = PlotlyJS.scatter(x=E, y=F, z=G, type="scatter3d", mode="lines", line=attr(color="blue", width=2))#Normal
    plot5 = PlotlyJS.scatter(x=EE, y=FF, z=GG, type="scatter3d", mode="lines", line=attr(color="blue", width=2))#Normal
    plot6 = PlotlyJS.scatter(x=[Cent[1]], y=[Cent[2]], z=[Cent[3]], type="scatter3d", mode="markers", marker=attr(size=2, color="blue", opacity=1))#Spindle Midpoint

    pp = PlotlyJS.plot([plot1, plot2, plot3], layout3_5)# plot both the surface and the centrosome
    resizefigtime = 0.75
    camera = attr(eye=attr(x=-1 * resizefigtime, y=-0.25 * resizefigtime, z=6.5 * resizefigtime))

    PlotlyJS.relayout!(pp, scene_camera=camera)

    display(pp)
end
