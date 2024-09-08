
### Create a files liost of all files in the folder and all its subfolders
function ListFiles(path::String)
    filelist = Vector{String}(undef, 0)
    for (root, dirs, files) in walkdir(path)
        for file in files
            push!(filelist, joinpath(root, file))
        end
    end
    return filelist
end

### Normalize angle so they are between 0 and 90 degrees
function CorrectAngles!(A)
    if A > 90
        A = 180 - A
    end
    return A
end

### Calculate angle in 3D between two 3D vectors 
function CalculerAngle(A, B)
    C = (acos(dot(A, B) / (norm(A) * norm(B)))) * 180 / pi
    C = CorrectAngles!(C)
    return C
end

### Calculate the Euclidean distance between two 3D popints
function Distance3D(A, B)
    distance = sqrt((A[1] - B[1])^2 + (A[2] - B[2])^2 + (A[3] - B[3])^2)
    return distance
end


### Remove the extension from the name of the file, eg "foo.csv" becomes "foo"
function RemoveExt(filename::String, ext::String)
    index = findlast(ext, filename)
    return filename[1:index[1]-1]
end

### 
function TxtfileList(filelist::Vector{String})
    return filter(x -> endswith(x, ".txt") && !(occursin("coords", x) || occursin("cellIDs", x) || occursin("Values", x) || occursin("console", x)), filelist)
end

function WrlfileList(filelist::Vector{String})
    return unique(filter(x -> endswith(x, ".wrl"), filelist))
end

function ZipfileList(filelist::Vector{String}, ends::String)
    return filter(x -> endswith(x, "$(ends).zip"), filelist)
end

function FetchDPaxis(filelist::Vector{String})
    DPaxis = Dict()
    for i in eachindex(filelist)
        if contains(filelist2[i], "_DPaxis.csv")
            movie = RemoveExt(basename(filelist[i]), "_DPaxis.csv")
            df = DataFrame(CSV.File(filelist[i]))
            select!(df, [:x1, :x2, :y1, :y2, :z1, :z2])
            df.X .= df[!, :x1] .- df[!, :x2]
            df.Y .= df[!, :y1] .- df[!, :y2]
            df.Z .= df[!, :z1] .- df[!, :z2]
            DPaxis[movie] = df
        end
    end
    return DPaxis
end

function FetchDTC(filelist::Vector{String})
    DTC = Dict()
    for i in eachindex(filelist)
        if contains(filelist[i], "_DTC.csv")
            movie = RemoveExt(basename(filelist[i]), "_DTC.csv")
            df = DataFrame(CSV.File(filelist[i]))
            select!(df, [:x1, :y1, :z1])
            DTC[movie] = df
        end
    end
    return DTC
end

function ReadTracks(txtfilepath::String)
    df = dropmissing(CSV.read(txtfilepath, DataFrame))
    filter!(:Label => !contains("ID"), df)
    filter!(:Label => contains("Cent"), df)
    return df
end

function ProcessLine1!(A)
    A = string(A)
    A = replace(A, "-1" => " ") # remove the signal between faces coordinates
    A = split(A, " ")
    A = chop.(A, head=0, tail=1)#remove comas
    A = filter!(x -> x !== nothing, tryparse.(Int64, A)) .+ 1 #Julia starts at 1, not 0
    A = collect(Iterators.partition(A, 3))
    return (A)
end

function ProcessLine2!(A)
    A = string(A)
    A = replace(A, "/" => " ")
    A = replace(A, "\\" => " ")
    A = replace(A, "," => " ")
    A = split(A, " ")
    A = filter!(x -> x !== nothing, tryparse.(Float64, A))
    return (A)
end

function process_wrl_files(A, D)
    Node_Count = 0
    B = Array{Int64}(undef, 0, 2)
    C = Array{Int64}(undef, 0, 2)
    for (y, line) in enumerate(eachline(A))
        if contains(line, D)
            Node_Count += 1
            fazer = [Node_Count y]
            B = vcat(B, fazer)
        elseif contains(line, "]") && Node_Count >= 1 && y - B[end, 2] > 20
            if size(C) == (0, 2)
                fazer = [Node_Count y]
                C = vcat(C, fazer)
            elseif y - C[end, 2] > 20 && C[end, 1] != Node_Count
                fazer = [Node_Count y]
                C = vcat(C, fazer)
            end
        end
    end
    return (B, C)
end


######## Fetch DPaxis
function FetchDPaxis(filelist::Vector{String})
    DPaxis = Dict()
    for i in eachindex(filelist)
        if contains(filelist[i], "_DPaxis.csv")
            gonad = RemoveExt(basename(filelist[i]), "_DPaxis.csv")
            df = DataFrame(CSV.File(filelist[i]))
            select!(df, [:x1, :x2, :y1, :y2, :z1, :z2])
            df.X .= df[!, :x1] .- df[!, :x2]
            df.Y .= df[!, :y1] .- df[!, :y2]
            df.Z .= df[!, :z1] .- df[!, :z2]
            DPaxis[gonad] = df
        end
    end
    return DPaxis
end


#########filtering Cell Shape DF
function FilteringCellShape(CellShape_df_final::DataFrame)
    filtered_df = copy(CellShape_df_final)
    filtered_df = filter(row ->
            row[:MeanCellVolume] > 75
                && row[:MeanCellVolume] < 250
                && row[:CellVolumeSlope] < 2
                && row[:CellVolumeSlope] > -2
                && row[:Duration] > 10
                && row[:VarianceCellVolume] < 400
                && row[:MeanDistanceDTC] < 100
                && row[:MeanAngleMeanCtoC] < 20
                && row[:stdAngleMeanCtoC] < 20,
        filtered_df)
    dropmissing!(filtered_df)
    return filtered_df
end

function calculate_normals_and_centroids(A, B, C)
    temporary_normals = Dict{Int64,Vector{Vector{Float64}}}()
    temporary_centroids = Dict{Int64,Vector{Vector{Float64}}}()
    temporary_areas = Dict{Int64,Vector{Float64}}()
    for x in 1:C
        centroids = Vector{Vector{Float64}}()
        normals = Vector{Vector{Float64}}()
        areas = Vector{Float64}()
        for y in 1:length(B[x])
            V1 = A[x][B[x][y][1]]
            V2 = A[x][B[x][y][2]]
            V3 = A[x][B[x][y][3]]
            VA = V2 - V1
            VB = V3 - V1
            normal = cross(VA, VB) / 2
            push!(normals, normal)
            centroid = [mean([V1[1], V2[1], V3[1]]), mean([V1[2], V2[2], V3[2]]), mean([V1[3], V2[3], V3[3]])]
            push!(centroids, centroid)
            area = 0.5(norm(cross(VA, VB)))
            push!(areas, area)
        end
        temporary_normals[x] = normals
        temporary_centroids[x] = centroids
        temporary_areas[x] = areas
    end
    return (temporary_normals, temporary_centroids, temporary_areas)
end

function CreateCelloutput(filelist::Vector{String})
    txtfilelist = TxtfileList(filelist)
    Celloutput = Dict{String,DataFrame}()
    unique_tracks = Dict()
    for x in eachindex(txtfilelist)
        txtfile = txtfilelist[x]
        movie = RemoveExt(basename(txtfile), ".txt")
        df = ReadTracks(txtfile)
        unique_tracks[x] = sort(unique(df[!, :Label]))
        for y in 1:2:length(unique_tracks[x])
            cell = unique_tracks[x][y][1:end-1]
            cell2 = cell[6:end]
            Cellule = string(movie, " Cell_", "$cell2")
            cella = string(cell, "a")
            cellb = string(cell, "b")
            cent_a = filter(:Label => isequal(cella), df)
            cent_b = filter(:Label => isequal(cellb), df)
            Celloutput["$Cellule"] = innerjoin(cent_a, cent_b, on=:FRAME, makeunique=true)
            spindle_mid_point = Vector{Matrix{Float64}}()
            spindle_lengths = Vector{Float64}()
            for z in 1:length(Celloutput["$Cellule"][!, :FRAME])
                push!(spindle_mid_point, [((Celloutput["$Cellule"][z, :POSITION_X] + Celloutput["$Cellule"][z, :POSITION_X_1]) / 2) ((Celloutput["$Cellule"][z, :POSITION_Y] + Celloutput["$Cellule"][z, :POSITION_Y_1]) / 2) ((Celloutput["$Cellule"][z, :POSITION_Z] + Celloutput["$Cellule"][z, :POSITION_Z_1]) / 2)])
                Celloutput["$Cellule"][z, :FRAME] = Celloutput["$Cellule"][z, :FRAME] + 1
                cent_a = [Celloutput["$Cellule"][z, :POSITION_X] Celloutput["$Cellule"][z, :POSITION_Y] Celloutput["$Cellule"][z, :POSITION_Z]]
                cent_b = [Celloutput["$Cellule"][z, :POSITION_X_1] Celloutput["$Cellule"][z, :POSITION_Y_1] Celloutput["$Cellule"][z, :POSITION_Z_1]]
                spindle_length = Distance3D(cent_a, cent_b)
                push!(spindle_lengths, spindle_length)
            end
            Celloutput["$Cellule"][!, :Spindle_Midpoint] = spindle_mid_point
            Celloutput["$Cellule"][!, :Spindle_Length] = spindle_lengths
        end
    end
    return Celloutput
end

findnearest(A::AbstractArray, t) = findmin(abs.(A .- t))[2]

function ParseWRL(wrl::String)
    io = open(wrl[1], "r")
    CoordsIndicesStart, CoordsIndicesEnd = process_wrl_files(io, "coordIndex")
    io = open(wrl[1], "r")
    foo2 = readlines(io)
    io = open(wrl[1], "r")
    PointsStart, PointsEnd = process_wrl_files(io, "point")
    number_of_nodes = length(CoordsIndicesStart[:, 1])
    number_of_nodes2 = length(CoordsIndicesEnd[:, 1])
    number_of_nodes3 = length(PointsStart[:, 1])
    number_of_nodes4 = length(PointsEnd[:, 1])
    coordIndices_temporary = Dict{Int64,Vector{SubArray{Int64,1,Vector{Int64},Tuple{UnitRange{Int64}},true}}}()
    points_temporary = Dict{Int64,Vector{SubArray{Union{Nothing,Float64},1,Vector{Union{Nothing,Float64}},Tuple{UnitRange{Int64}},true}}}()
    Normals_temporary = Dict{Int64,Vector{Vector{Float64}}}()
    Areas_temporary = Dict{Int64,Vector{Float64}}()
    Centroids_temporary = Dict{Int64,Vector{Vector{Float64}}}()
    if number_of_nodes == number_of_nodes2 && number_of_nodes3 == number_of_nodes4 && number_of_nodes == number_of_nodes3
        for y in 1:number_of_nodes
            coordIndices_node = foo2[CoordsIndicesStart[y, 2]:CoordsIndicesEnd[y, 2]]
            coordIndices_node = process_line_i!(coordIndices_node)
            coordIndices_temporary[y] = coordIndices_node
            points_node = foo2[PointsStart[y, 2]:PointsEnd[y, 2]]
            points_node = process_line_ii!(points_node)
            points_node = collect(Iterators.partition(points_node, 3))
            points_temporary[y] = points_node
        end
        Normals_temporary, Centroids_temporary, Areas_temporary = calculate_normals_and_centroids(points_temporary, coordIndices_temporary, number_of_nodes)
    end
    close(io)
    return number_of_nodes, number_of_nodes2, number_of_nodes3, number_of_nodes4, Normals_temporary, Centroids_temporary, Areas_temporary
end

function CelloutputDict2DF(Celloutput)
    Celloutput = sort(Celloutput)
    k = sort(collect(keys(Celloutput)))
    Celloutputdf = DataFrame()
    for i in eachindex(k)
        df = Celloutput[k[i]]
        IndexCell = findfirst(" Cell", k[i])
        CelluleNumber = k[i][IndexCell[1]+1:end]
        println(CelluleNumber)
        GonadName = k[i][1:IndexCell[1]]
        cell = fill(CelluleNumber, (size(Celloutput[k[i]])[1]))
        gonad = fill(GonadName, (size(Celloutput[k[i]])[1]))
        df[!, :cell] = cell
        df[!, :gonad] = gonad
        if hasproperty(df, :AngleSpindleRachis) == false
            AngleSpindleRachis = fill(missing, (size(Celloutput[k[i]])[1]))
            df[!, :AngleSpindleRachis] = AngleSpindleRachis
        end
        if hasproperty(df, :XNormal) == false
            XNormal = fill(missing, (size(Celloutput[k[i]])[1]))
            YNormal = fill(missing, (size(Celloutput[k[i]])[1]))
            ZNormal = fill(missing, (size(Celloutput[k[i]])[1]))
            df[!, :XNormal] = XNormal
            df[!, :YNormal] = YNormal
            df[!, :ZNormal] = ZNormal
        end
        Celloutputdf = vcat(Celloutputdf, df)
    end
    Celloutputdf = DataFramesMeta.groupby(Celloutputdf, [:gonad, :cell])
 
    Celloutputdf = @combine(Celloutputdf, 
    :gonad = mode(:gonad),
    :cell = mode(:cell),
    :Time = [:POSITION_T],
    :Frame = [:FRAME],
    :Duration = length(:FRAME),
    :CentaX = [:POSITION_X],
    :CentaY = [:POSITION_Y],
    :CentaZ = [:POSITION_Z],
    :CentbX = [:POSITION_X_1],
    :CentbY = [:POSITION_Y_1],
    :CentbZ = [:POSITION_Z_1],
    :Spindle_Midpoint = [:Spindle_Midpoint],
    :Spindle_Length = [:Spindle_Length],
    :XNormal = [:XNormal], 
    :YNormal = [:YNormal],
    :ZNormal = [:ZNormal],
    :AngleSpindleRachis = [:AngleSpindleRachis],
    )
    return Celloutputdf
    #CSV.write("Celloutput.csv", Celloutputdf)
end

function AnglesMeanCtoC(df::DataFrame)
    AngleMeanCtoC = Vector{Vector{Union{Float64,Missing}}}(undef, nrow(df))
    MeanAngleMeanCtoC = Vector{Union{Float64,Missing}}(undef, nrow(df))
    stdAngleMeanCtoC = Vector{Union{Float64,Missing}}(undef, nrow(df))
    for (i, row) in enumerate(eachrow(df))
        MeanC = [row[:MeanCellEllipsoidAxisCX], row[:MeanCellEllipsoidAxisCY], row[:MeanCellEllipsoidAxisCZ]]
        AngleMeanCtoCT = Vector{Union{Float64}}(undef, length(row[:CellEllipsoidAxisCX]))
        for j in 1:length(row[:CellEllipsoidAxisCX])
            C = [row[:CellEllipsoidAxisCX][j], row[:CellEllipsoidAxisCY][j], row[:CellEllipsoidAxisCZ][j]]
            AngleMeanCtoCT[j] = CalculerAngle(MeanC, C)
        end
        AngleMeanCtoC[i] = AngleMeanCtoCT
        MeanAngleMeanCtoC[i] = mean(AngleMeanCtoCT)
        stdAngleMeanCtoC[i] = std(AngleMeanCtoCT)
    end
    df[!, :AngleMeanCtoC] = AngleMeanCtoC
    df[!, :MeanAngleMeanCtoC] = MeanAngleMeanCtoC
    df[!, :stdAngleMeanCtoC] = stdAngleMeanCtoC
    return df
end

function ProjectVectorSurface(v::Vector{Float64}, n::Vector{Float64})
    # Normalize the normal vector
    n_normalized = normalize(n)

    # Compute the projection of v onto n
    projection_onto_n = dot(v, n_normalized) * n_normalized

    # Subtract the projection from the original vector to get the component in the plane
    projection_onto_plane = v - projection_onto_n

    return projection_onto_plane
end

function FractionCumSum(fazer::Vector)
    cumsum(sort(fazer)) ./ sum((fazer))
end
#########filtering Cell Shape DF
function FilteringCelloutput(CellShape_df_final::DataFrame)
    filtered_df = copy(CellShape_df_final)
    filtered_df = filter(row ->
            row[:MeanCellVolume] > 75
                && row[:MeanCellVolume] < 250
                && row[:CellVolumeSlope] < 2
                && row[:CellVolumeSlope] > -2
                && row[:VarianceCellVolume] < 400
                && row[:MeanDistanceDTC] < 100
                && row[:MeanAngleMeanCtoC] < 20
                && row[:stdAngleMeanCtoC] < 20,
        filtered_df)
    return filtered_df
end

function KeepBigCells(CellShapeCleaned::DataFrame, Celloutputiii::DataFrame)
    filtered_df = copy(CellShapeCleaned)
    MeanCellVolumeMitosis = mean(Celloutputiii[!, :MeanCellVolume])
    StdCellVolumeMitosis = std(Celloutputiii[!, :MeanCellVolume])
    filtered_df = filter(row ->
            row[:MeanCellVolume] > MeanCellVolumeMitosis - (2 * StdCellVolumeMitosis)
            &&
                row[:MeanCellVolume] < MeanCellVolumeMitosis + (2 * StdCellVolumeMitosis),
        filtered_df)
    return filtered_df
end

function ProjectVonU(V, U)
    û = U / norm(U)
    pV = V .- (dot(û, V) * û)
    return pV
end

function PolarHist(angles, figure, title, colour, transparancy, graphposition)
    anglesRad = deg2rad.(angles)
    binsRad = deg2rad.(collect(0:15:90))
    h = fit(Histogram, anglesRad, binsRad)
    x = collect(h.edges[1])
    y = (h.weights ./ sum(h.weights)) .* 100
    ax = PolarAxis(figure[1, graphposition],
        title="$title",
        thetaticks=((collect(0:15:90)) .* pi ./ 180, string.(string.(collect(0:15:90)), "°")),
        rticks=collect(0:10:ceil(maximum(y) / 10)*10)
    )
    thetalims!(ax, 0, pi / 2)
    rlims!(ax, 0, ceil(maximum(y) / 10) * 10)
    for i in eachindex(y)
        if y[i] > 0
            CairoMakie.poly!(ax, [(0.0, 0.0), (x[i], y[i]), (x[i+1], y[i])], [1, 2, 3], color=Makie.wong_colors()[colour], alpha=transparancy, strokewidth=2, strokecolor=:black)
        end
    end
end

function ParseDataframe(df)
    namesdf = names(df)
    for name in namesdf
        if sum(typeof(df[!, name]) .== [Vector{Union{Missing,Int64}}, Vector{Int64}, Vector{Union{Missing,Float64}}, Vector{Float64}]) == 0
            println(name)
            println(typeof(df[!, name]))
            fazel = []
            fazel1 = Vector{Vector{Union{Missing,Float64}}}(undef, 0)
            fazel2 = Vector{Vector{Union{Missing,Int64}}}(undef, 0)
            fazel3 = []

            for i in 1:nrow(df)
                if startswith(string(df[i, name]), "Union{Missing, Float64}")
                    push!(fazel1, passmissing(parse).(Float64, replace(strip.(split(df[i, name][25:end-1], ",")), "missing" => missing)))
                elseif startswith(string(df[i, name]), "Union{Missing, Int64}")
                    push!(fazel2, passmissing(parse).(Int64, replace(strip.(split(df[i, name][23:end-1], ",")), "missing" => missing)))
                elseif startswith(string(df[i, name]), "[")
                    push!(fazel1, passmissing(parse).(Float64, replace(strip.(split(df[i, name][2:end-1], ",")), "missing" => missing)))
                else
                    push!(fazel3, df[i, name])
                end
            end
            SizeFazel = Vector{Int64}(undef, 3)
            SizeFazel[1] = length(fazel1)
            SizeFazel[2] = length(fazel2)
            SizeFazel[3] = length(fazel3)
            IndexSF = findmax(SizeFazel)
            if IndexSF[2] == 1
                fazel = fazel1
            elseif IndexSF[2] == 2
                fazel = fazel2
            else
                fazel = fazel3
            end
            df[!, name] = fazel
            fazel = []
            fazel1 = []
            fazel2 = []
        end
    end
    return df
end

function RemoveMitoticTrackID(Celloutput, CellShape)
    df = copy(CellShape)
    for i in 1:nrow(Celloutput)
        gonad = strip(Celloutput[i, :gonad])
        for j in 1:Celloutput[i, :Duration]
            TrackID = Celloutput[i, :CellShapeTrackID][j]
            if ismissing(TrackID) == false
                df = filter(row -> !(row.TrackID == TrackID && strip(row.Gonad) == gonad), df)
            end
        end
    end
    return df
end

function CreateInterphaseCellsDF(filtered_df)
    InterphaseCellsDF = DataFrame()
    Gonad = []
    TrackID = []
    Type = []
    Time = []
    Duration = []
    CellPositionX = []
    CellPositionY = []
    CellPositionZ = []
    CellEllipsoidAxisCX = []
    CellEllipsoidAxisCY = []
    CellEllipsoidAxisCZ = []
    MeanDistanceDTC = []
    MeanCellVolume = []
    MeanOrientationCtoDP = []
    MeanProlateIndex = []

    for i in 1:nrow(filtered_df)
        Gonad = push!(Gonad, filtered_df[i, :Gonad])
        TrackID = push!(TrackID, filtered_df[i, :TrackID])
        Type = push!(Type, filtered_df[i, :Type])
        Time = push!(Time, filtered_df[i, :Time][1])
        Duration = push!(Duration, filtered_df[i, :Duration])
        CellPositionX = push!(CellPositionX, filtered_df[i, :CellPositionX][1])
        CellPositionY = push!(CellPositionY, filtered_df[i, :CellPositionY][1])
        CellPositionZ = push!(CellPositionZ, filtered_df[i, :CellPositionZ][1])
        CellEllipsoidAxisCX = push!(CellEllipsoidAxisCX, filtered_df[i, :CellEllipsoidAxisCX][1])
        CellEllipsoidAxisCY = push!(CellEllipsoidAxisCY, filtered_df[i, :CellEllipsoidAxisCY][1])
        CellEllipsoidAxisCZ = push!(CellEllipsoidAxisCZ, filtered_df[i, :CellEllipsoidAxisCZ][1])
        MeanDistanceDTC = push!(MeanDistanceDTC, filtered_df[i, :MeanDistanceDTC])
        MeanCellVolume = push!(MeanCellVolume, filtered_df[i, :MeanCellVolume])
        MeanOrientationCtoDP = push!(MeanOrientationCtoDP, filtered_df[i, :MeanOrientationCtoDP])
        MeanProlateIndex = push!(MeanProlateIndex, filtered_df[i, :MeanProlateIndex])
    end

    InterphaseCellsDF[!, :Gonad] = Gonad
    InterphaseCellsDF[!, :TrackID] = TrackID
    InterphaseCellsDF[!, :Type] = Type
    InterphaseCellsDF[!, :Time] = Time
    InterphaseCellsDF[!, :Duration] = Duration
    InterphaseCellsDF[!, :CellPositionX] = CellPositionX
    InterphaseCellsDF[!, :CellPositionY] = CellPositionY
    InterphaseCellsDF[!, :CellPositionZ] = CellPositionZ
    InterphaseCellsDF[!, :CellEllipsoidAxisCX] = CellEllipsoidAxisCX
    InterphaseCellsDF[!, :CellEllipsoidAxisCY] = CellEllipsoidAxisCY
    InterphaseCellsDF[!, :CellEllipsoidAxisCZ] = CellEllipsoidAxisCZ
    InterphaseCellsDF[!, :MeanDistanceDTC] = MeanDistanceDTC
    InterphaseCellsDF[!, :MeanCellVolume] = MeanCellVolume
    InterphaseCellsDF[!, :MeanOrientationCtoDP] = MeanOrientationCtoDP
    InterphaseCellsDF[!, :MeanProlateIndex] = MeanProlateIndex

    return InterphaseCellsDF
end

function AddDPaxisData(df)
    df[!, :DPaxisX] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    df[!, :DPaxisY] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    df[!, :DPaxisZ] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    df[!, :AngleSpindleDPaxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(df))
    df[!, :MeanAngleSpindleDPaxis] = Vector{Union{Missing,Float64}}(undef, nrow(df))

    for i in 1:nrow(df)
        println(i)
        df[i, :AngleSpindleDPaxis] = Vector{Union{Missing,Float64}}(undef, df[i, :Duration])
        gonad = strip.(df[i, :gonad])
        DPaxisVect = [missing, missing, missing]
        if haskey(DPaxis, gonad) == true
            DPaxisVect = [DPaxis[gonad][!, :X][1], DPaxis[gonad][!, :Y][1], DPaxis[gonad][!, :Z][1]]
            df[i, :DPaxisX] = DPaxisVect[1]
            df[i, :DPaxisY] = DPaxisVect[2]
            df[i, :DPaxisZ] = DPaxisVect[3]
        end
        for j in 1:df[i, :Duration]
            Spindle = [df[i, :CentaX][j] - df[i, :CentbX][j], df[i, :CentaY][j] - df[i, :CentbY][j], df[i, :CentaZ][j] - df[i, :CentbZ][j]]
            df[i, :AngleSpindleDPaxis][j] = CalculerAngle(DPaxisVect, Spindle)
        end
        IndexAnaphaseonset = findfirst(isequal(df[i, :AnaphaseOnset]), df[i, :Frame])
        if length(df[i, :Frame][1:IndexAnaphaseonset]) >= 5 && IndexAnaphaseonset + 1 <= length(df[i, :Frame])
            df[i, :MeanAngleSpindleDPaxis] = mean(df[i, :AngleSpindleDPaxis][IndexAnaphaseonset-4:IndexAnaphaseonset])
        end
    end
    return df
end

function AddDPaxisData(df)
    df[!, :DPaxisX] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    df[!, :DPaxisY] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    df[!, :DPaxisZ] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    df[!, :AngleSpindleDPaxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(df))
    df[!, :AngleLongAxisDPaxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(df))
    df[!, :MeanAngleSpindleDPaxis] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    df[!, :MeanAngleLongAxisDPaxis] = Vector{Union{Missing,Float64}}(undef, nrow(df))

    for i in 1:nrow(df)
        println(i)
        df[i, :AngleSpindleDPaxis] = Vector{Union{Missing,Float64}}(undef, df[i, :Duration])
        df[i, :AngleLongAxisDPaxis] = Vector{Union{Missing,Float64}}(undef, df[i, :Duration])
        gonad = strip.(df[i, :gonad])
        DPaxisVect = [missing, missing, missing]
        if haskey(DPaxis, gonad) == true
            DPaxisVect = [DPaxis[gonad][!, :X][1], DPaxis[gonad][!, :Y][1], DPaxis[gonad][!, :Z][1]]
            df[i, :DPaxisX] = DPaxisVect[1]
            df[i, :DPaxisY] = DPaxisVect[2]
            df[i, :DPaxisZ] = DPaxisVect[3]
        end
        for j in 1:df[i, :Duration]
            Spindle = [df[i, :CentaX][j] - df[i, :CentbX][j], df[i, :CentaY][j] - df[i, :CentbY][j], df[i, :CentaZ][j] - df[i, :CentbZ][j]]
            LongAxis = [df[i, :CellEllipsoidAxisCX][j], df[i, :CellEllipsoidAxisCY][j], df[i, :CellEllipsoidAxisCZ][j]]
            df[i, :AngleSpindleDPaxis][j] = CalculerAngle(DPaxisVect, Spindle)
            if sum(ismissing.(LongAxis)) == 0
                df[i, :AngleLongAxisDPaxis][j] = CalculerAngle(DPaxisVect, LongAxis)
            end
        end
        IndexAnaphaseonset = findfirst(isequal(df[i, :AnaphaseOnset]), df[i, :Frame])
        if length(df[i, :Frame][1:IndexAnaphaseonset]) >= 5 && IndexAnaphaseonset + 1 <= length(df[i, :Frame])
            df[i, :MeanAngleSpindleDPaxis] = mean(df[i, :AngleSpindleDPaxis][IndexAnaphaseonset-4:IndexAnaphaseonset])
            df[i, :MeanAngleLongAxisDPaxis] = mean(df[i, :AngleLongAxisDPaxis][IndexAnaphaseonset-4:IndexAnaphaseonset])
        end
    end
    return df
end

function AddDPaxisDataInterphase(df)
    df[!, :DPaxisX] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    df[!, :DPaxisY] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    df[!, :DPaxisZ] = Vector{Union{Missing,Float64}}(undef, nrow(df))

    df[!, :AngleLongAxisDPaxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(df))

    df[!, :MeanAngleLongAxisDPaxis] = Vector{Union{Missing,Float64}}(undef, nrow(df))

    for i in 1:nrow(df)
        println(i)

        df[i, :AngleLongAxisDPaxis] = Vector{Union{Missing,Float64}}(undef, df[i, :Duration])
        gonad = strip.(df[i, :Gonad])
        DPaxisVect = [missing, missing, missing]
        if haskey(DPaxis, gonad) == true
            DPaxisVect = [DPaxis[gonad][!, :X][1], DPaxis[gonad][!, :Y][1], DPaxis[gonad][!, :Z][1]]
            df[i, :DPaxisX] = DPaxisVect[1]
            df[i, :DPaxisY] = DPaxisVect[2]
            df[i, :DPaxisZ] = DPaxisVect[3]
        end
        for j in 1:df[i, :Duration]

            LongAxis = [df[i, :CellEllipsoidAxisCX][j], df[i, :CellEllipsoidAxisCY][j], df[i, :CellEllipsoidAxisCZ][j]]

            if sum(ismissing.(LongAxis)) == 0
                df[i, :AngleLongAxisDPaxis][j] = CalculerAngle(DPaxisVect, LongAxis)
            end
        end
        df[i, :MeanAngleLongAxisDPaxis] = mean(df[i, :AngleLongAxisDPaxis][1:10])
    end
    return df
end

function SaveInterphaseDfPublication(df, path, name)
    dfp = select(df,
        [
            :Gonad,
            :TrackID,
            :MeanDistanceDTC,
            :MeanCellVolume,
            :MeanProlateIndex,
            :MeanAnglesLongAxisRachis,
            :MeanAngleLongAxisDPaxis, 
            :MeanAngleProjectLongestAxisDP,
        ]
    )
    dfp.Type = fill("Interphase", nrow(dfp))
    rename!(dfp, Dict( 
        :MeanDistanceDTC => "Distance_DTC",
        :MeanCellVolume => "Cell_volume",
        :MeanProlateIndex => "Prolate_index",
        :MeanAnglesLongAxisRachis => "Angle_long_axis_rachis",
        :MeanAngleLongAxisDPaxis => "Angle_long_axis_DP_axis",
        :MeanAngleProjectLongestAxisDP => "Angle_longest_projected_axis_DP_axis"
    ))
    CSV.write(joinpath(path, "$(name)CellShapeData.csv"), dfp)
end


function SaveMitosisDfPublication(df)
    dfp = select(df,
        [
            :gonad,
            :cell,
            :TrackID,
            :MeanDistanceDTC,
            :MeanCellVolume,
            :MeanProlateIndex,
            :MeanAngleLongAxisRachis,
            :MeanAngleLongAxisDPaxis,
            :MeanAngleSpindleRachisAnaphaseOnset,
            :MeanAngleSpindleLongAxisAnaphaseOnset,
            :MeanAngleSpindleShortAxis,
            :MeanAngleSpindleDPaxis,
            :MeanAngleProjectSpindleProjectLongAxis,
            :MeanAngleProjectSpindleProjectDPaxis,
            :MeanAngleProjectLongestAxisDP,
            ]
    )
    dfp.Type = fill("Mitosis", nrow(dfp))
    rename!(dfp, Dict(
        :gonad => "Gonad",
        :cell => "Cell", 
        :MeanDistanceDTC => "Distance_DTC",
        :MeanCellVolume => "Cell_volume",
        :MeanProlateIndex => "Prolate_index",
        :MeanAngleLongAxisRachis => "Angle_long_axis_rachis",
        :MeanAngleLongAxisDPaxis => "Angle_long_axi_DP_axis",
        :MeanAngleSpindleRachisAnaphaseOnset => "Angle_spindle_rachis",
        :MeanAngleSpindleLongAxisAnaphaseOnset => "Angle_spindle_long_axis",
        :MeanAngleSpindleShortAxis => "Angle_spindle_short_axis",
        :MeanAngleSpindleDPaxis => "Angle_spindle_DP_axis",
        :MeanAngleProjectSpindleProjectLongAxis => "Angle_project_spindle_project_long_axis",
        :MeanAngleProjectSpindleProjectDPaxis => "Angle_project_spindle_DP_axis",
        :MeanAngleProjectLongestAxisDP => "Angle_project_long_axis_DP_axis"
    ))
    #CSV.write(joinpath(path, "MitosisCellShapeData.csv"), dfp)
    return dfp
end

function AddSpindleSmallAxis(df)
    df[!, :AngleSpindleShortAxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(df))
    df[!, :MeanAngleSpindleShortAxis] = Vector{Union{Missing,Float64}}(undef, nrow(df))
    for i in 1:nrow(df)
        if ismissing(df[i, :AnaphaseOnset]) == false
            IndexAnaphaseonset = findfirst(isequal(df[i, :AnaphaseOnset]), df[i, :Frame])
            if length(df[i, :Frame][1:IndexAnaphaseonset]) >= 5 && IndexAnaphaseonset + 1 <= length(df[i, :Frame])
                df[i, :AngleSpindleShortAxis] = Vector{Union{Missing,Float64}}(undef, df[i, :Duration])
                for j in 1:df[i, :Duration]
                    ShortAxis = [df[i, :CellEllipsoidAxisAX][j], df[i, :CellEllipsoidAxisAY][j], df[i, :CellEllipsoidAxisAZ][j]] .* df[i, :CellEllipsoidAxisLengthA][j]
                    Spindle = [df[i, :CentaX][j] - df[i, :CentbX][j], df[i, :CentaY][j] - df[i, :CentbY][j], df[i, :CentaZ][j] - df[i, :CentbZ][j]]
                    if sum(ismissing.(Spindle)) == 0 && sum(ismissing.(ShortAxis)) == 0
                        df[i, :AngleSpindleShortAxis][j] = CalculerAngle(ShortAxis, Spindle)
                    end
                end
                df[i, :MeanAngleSpindleShortAxis] = mean(df[i, :AngleSpindleShortAxis][IndexAnaphaseonset-4:IndexAnaphaseonset])
            end
        end
    end
    return df
end