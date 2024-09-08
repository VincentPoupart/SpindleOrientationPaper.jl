######## Fetch Data from Imaris

function FetchImarisStats(filelist::Vector{String}, DTC)
    Cell_Ellipsoid_Axis_C_df = DataFrame()
    Cell_Ellipsoid_Axis_B_df = DataFrame()
    Cell_Ellipsoid_Axis_A_df = DataFrame()
    Cell_Ellipsoid_Axis_Length_C_df = DataFrame()
    Cell_Ellipsoid_Axis_Length_B_df = DataFrame()
    Cell_Ellipsoid_Axis_Length_A_df = DataFrame()
    Cell_Volume_df = DataFrame()
    Cell_Position_df = DataFrame()
    Cell_Prolate_df = DataFrame()
    for i in eachindex(filelist)
        if contains(filelist[i], "Track") == false && contains(filelist[i], "Overall") == false
            if contains(filelist[i], "Cell_Ellipsoid_Axis_C")
                df = DataFrame(CSV.File(filelist[i], header=["CellEllipsoidAxisCX", "CellEllipsoidAxisCY", "CellEllipsoidAxisCZ", "Unit", "Category", "Collection", "Time", "TrackID", "ID", "Column10"], skipto=5, delim=','))
                gonad = RemoveExt(basename(filelist[i]), "_Cell_Ellipsoid_Axis_C.csv")
                g = fill(gonad, nrow(df))
                h = fill("Cell", nrow(df))
                df[!, :Gonad] = g
                df[!, :Type] = h
                NamesDF = names(df)
                IndexColRemove = findall(contains("Column"), NamesDF)
                select!(df, Not(NamesDF[IndexColRemove]))
                select!(df, Not([:ID, :Unit, :Category, :Collection]))
                dropmissing!(df)
                angles = Vector{Union{Float64,Missing}}(undef, nrow(df))
                if haskey(DPaxis, gonad) == true
                    DPaxisVect = [DPaxis[gonad][!, :X], DPaxis[gonad][!, :Y], DPaxis[gonad][!, :Z]]
                    for j in 1:nrow(df)
                        AxisCVect = [df[j, :CellEllipsoidAxisCX], df[j, :CellEllipsoidAxisCY], df[j, :CellEllipsoidAxisCZ]]
                        angles[j] = CalculerAngle(DPaxisVect, AxisCVect)

                    end
                end
                df[!, :OrientationCtoDP] = angles
                Cell_Ellipsoid_Axis_C_df = vcat(Cell_Ellipsoid_Axis_C_df, df)
                df = []
            elseif contains(filelist[i], "Cell_Ellipsoid_Axis_B")
                df = DataFrame(CSV.File(filelist[i], header=["CellEllipsoidAxisBX", "CellEllipsoidAxisBY", "CellEllipsoidAxisBZ", "Unit", "Category", "Collection", "Time", "TrackID", "ID", "Column10"], skipto=5, delim=','))
                gonad = RemoveExt(basename(filelist[i]), "_Cell_Ellipsoid_Axis_B.csv")
                g = fill(gonad, nrow(df))
                h = fill("Cell", nrow(df))
                df[!, :Gonad] = g
                df[!, :Type] = h
                NamesDF = names(df)
                IndexColRemove = findall(contains("Column"), NamesDF)
                select!(df, Not(NamesDF[IndexColRemove]))
                select!(df, Not([:ID, :Unit, :Category, :Collection]))
                dropmissing!(df)
                angles = Vector{Union{Float64,Missing}}(undef, nrow(df))
                if haskey(DPaxis, gonad) == true
                    DPaxisVect = [DPaxis[gonad][!, :X], DPaxis[gonad][!, :Y], DPaxis[gonad][!, :Z]]
                    for j in 1:nrow(df)
                        AxisBVect = [df[j, :CellEllipsoidAxisBX], df[j, :CellEllipsoidAxisBY], df[j, :CellEllipsoidAxisBZ]]
                        angles[j] = CalculerAngle(DPaxisVect, AxisBVect)

                    end
                end
                df[!, :OrientationBtoDP] = angles
                Cell_Ellipsoid_Axis_B_df = vcat(Cell_Ellipsoid_Axis_B_df, df)
                df = []
            elseif contains(filelist[i], "Cell_Ellipsoid_Axis_A")
                df = DataFrame(CSV.File(filelist[i], header=["CellEllipsoidAxisAX", "CellEllipsoidAxisAY", "CellEllipsoidAxisAZ", "Unit", "Category", "Collection", "Time", "TrackID", "ID", "Column10"], skipto=5, delim=','))
                gonad = RemoveExt(basename(filelist[i]), "_Cell_Ellipsoid_Axis_A.csv")
                g = fill(gonad, nrow(df))
                h = fill("Cell", nrow(df))
                df[!, :Gonad] = g
                df[!, :Type] = h
                NamesDF = names(df)
                IndexColRemove = findall(contains("Column"), NamesDF)
                select!(df, Not(NamesDF[IndexColRemove]))
                select!(df, Not([:ID, :Unit, :Category, :Collection]))
                dropmissing!(df)
                angles = Vector{Union{Float64,Missing}}(undef, nrow(df))
                if haskey(DPaxis, gonad) == true
                    DPaxisVect = [DPaxis[gonad][!, :X], DPaxis[gonad][!, :Y], DPaxis[gonad][!, :Z]]
                    for j in 1:nrow(df)
                        AxisAVect = [df[j, :CellEllipsoidAxisAX], df[j, :CellEllipsoidAxisAY], df[j, :CellEllipsoidAxisAZ]]
                        angles[j] = CalculerAngle(DPaxisVect, AxisAVect)

                    end
                end
                df[!, :OrientationAtoDP] = angles
                Cell_Ellipsoid_Axis_A_df = vcat(Cell_Ellipsoid_Axis_A_df, df)
                df = []
            elseif contains(filelist[i], "Cell_Ellipsoid_Axis_Length_C")
                df = DataFrame(CSV.File(filelist[i], header=["CellEllipsoidAxisLengthC", "Unit", "Category", "Time", "TrackID", "ID", "Column7"], skipto=5, delim=','))
                gonad = RemoveExt(basename(filelist[i]), "_Cell_Ellipsoid_Axis_Length_C.csv")
                g = fill(gonad, nrow(df))
                h = fill("Cell", nrow(df))
                df[!, :Gonad] = g
                df[!, :Type] = h
                NamesDF = names(df)
                IndexColRemove = findall(contains("Column"), NamesDF)
                select!(df, Not(NamesDF[IndexColRemove]))
                select!(df, Not([:ID, :Unit, :Category]))
                dropmissing!(df)
                Cell_Ellipsoid_Axis_Length_C_df = vcat(Cell_Ellipsoid_Axis_Length_C_df, df)
                df = []
            elseif contains(filelist[i], "Cell_Ellipsoid_Axis_Length_B")
                df = DataFrame(CSV.File(filelist[i], header=["CellEllipsoidAxisLengthB", "Unit", "Category", "Time", "TrackID", "ID", "Column7"], skipto=5, delim=','))
                gonad = RemoveExt(basename(filelist[i]), "_Cell_Ellipsoid_Axis_Length_B.csv")
                g = fill(gonad, nrow(df))
                h = fill("Cell", nrow(df))
                df[!, :Gonad] = g
                df[!, :Type] = h
                NamesDF = names(df)
                IndexColRemove = findall(contains("Column"), NamesDF)
                select!(df, Not(NamesDF[IndexColRemove]))
                select!(df, Not([:ID, :Unit, :Category]))
                dropmissing!(df)
                Cell_Ellipsoid_Axis_Length_B_df = vcat(Cell_Ellipsoid_Axis_Length_B_df, df)
                df = []
            elseif contains(filelist[i], "Cell_Ellipsoid_Axis_Length_A")
                df = DataFrame(CSV.File(filelist[i], header=["CellEllipsoidAxisLengthA", "Unit", "Category", "Time", "TrackID", "ID", "Column7"], skipto=5, delim=','))
                gonad = RemoveExt(basename(filelist[i]), "_Cell_Ellipsoid_Axis_Length_A.csv")
                g = fill(gonad, nrow(df))
                h = fill("Cell", nrow(df))
                df[!, :Gonad] = g
                df[!, :Type] = h
                NamesDF = names(df)
                IndexColRemove = findall(contains("Column"), NamesDF)
                select!(df, Not(NamesDF[IndexColRemove]))
                select!(df, Not([:ID, :Unit, :Category]))
                dropmissing!(df)
                Cell_Ellipsoid_Axis_Length_A_df = vcat(Cell_Ellipsoid_Axis_Length_A_df, df)
                df = []
            elseif contains(filelist[i], "Cell_Volume")
                df = DataFrame(CSV.File(filelist[i], header=["CellVolume", "Unit", "Category", "Time", "TrackID", "ID", "Column7"], skipto=5, delim=','))
                gonad = RemoveExt(basename(filelist[i]), "_Cell_Volume.csv")
                g = fill(gonad, nrow(df))
                df[!, :Gonad] = g
                NamesDF = names(df)
                IndexColRemove = findall(contains("Column"), NamesDF)
                select!(df, Not(NamesDF[IndexColRemove]))
                select!(df, Not([:ID, :Unit, :Category]))
                dropmissing!(df)
                Cell_Volume_df = vcat(Cell_Volume_df, df)
                df = []
            elseif contains(filelist[i], "Cell_Ellipticity_(prolate)")
                df = DataFrame(CSV.File(filelist[i], header=["Cell_Ellipticity_Prolate", "Unit", "Category", "Time", "TrackID", "ID", "Column7"], skipto=5, delim=','))
                gonad = RemoveExt(basename(filelist[i]), "_Cell_Ellipticity_(prolate).csv")
                g = fill(gonad, nrow(df))
                df[!, :Gonad] = g
                NamesDF = names(df)
                IndexColRemove = findall(contains("Column"), NamesDF)
                select!(df, Not(NamesDF[IndexColRemove]))
                select!(df, Not([:ID, :Unit, :Category]))
                dropmissing!(df)
                Cell_Prolate_df = vcat(Cell_Prolate_df, df)
                df = []
            elseif contains(filelist[i], "Cell_Position")
                df = DataFrame(CSV.File(filelist[i], header=["CellPositionX", "CellPositionY", "CellPositionZ", "Unit", "Category", "Collection", "Time", "TrackID", "ID", "Column10"], skipto=5))
                gonad = RemoveExt(basename(filelist[i]), "_Cell_Position.csv")
                g = fill(gonad, nrow(df))
                df[!, :Gonad] = g
                NamesDF = names(df)
                IndexColRemove = findall(contains("Column"), NamesDF)
                select!(df, Not(NamesDF[IndexColRemove]))
                select!(df, Not([:ID, :Unit, :Category, :Collection]))
                dropmissing!(df)
                distanceDTC = Vector{Union{Float64,Missing}}(undef, nrow(df))
                if haskey(DTC, gonad) == true
                    DTC_Position = [DTC[gonad][1, :x1], DTC[gonad][1, :y1], DTC[gonad][1, :z1]]
                    for j in 1:nrow(df)
                        Cell_Position = [df[j, :CellPositionX], df[j, :CellPositionY], df[j, :CellPositionZ]]
                        distanceDTC[j] = Distance3D(DTC_Position, Cell_Position)
                    end
                end
                df[!, :DistanceDTC] = distanceDTC
                Cell_Position_df = vcat(Cell_Position_df, df)
                df = []
            end
        end
    end
    Cell_Shape_df = leftjoin(Cell_Ellipsoid_Axis_A_df, Cell_Ellipsoid_Axis_B_df; on=[:Gonad, :Time, :TrackID, :Type])
    Cell_Shape_df = leftjoin(Cell_Shape_df, Cell_Ellipsoid_Axis_C_df; on=[:Gonad, :Time, :TrackID, :Type])
    Cell_Shape_df = leftjoin(Cell_Shape_df, Cell_Ellipsoid_Axis_Length_C_df; on=[:Gonad, :Time, :TrackID, :Type])
    Cell_Shape_df = leftjoin(Cell_Shape_df, Cell_Ellipsoid_Axis_Length_B_df; on=[:Gonad, :Time, :TrackID, :Type])
    Cell_Shape_df = leftjoin(Cell_Shape_df, Cell_Ellipsoid_Axis_Length_A_df; on=[:Gonad, :Time, :TrackID, :Type])
    Cell_Shape_df = leftjoin(Cell_Shape_df, Cell_Volume_df; on=[:Gonad, :Time, :TrackID])
    Cell_Shape_df = leftjoin(Cell_Shape_df, Cell_Position_df; on=[:Gonad, :Time, :TrackID])
    Cell_Shape_df = leftjoin(Cell_Shape_df, Cell_Prolate_df; on=[:Gonad, :Time, :TrackID])
    dropmissing!(Cell_Shape_df)
    Cell_Shape_groupeddf = DataFramesMeta.groupby(Cell_Shape_df, [:Gonad, :TrackID])

    CellShape_df_final = @combine(Cell_Shape_groupeddf,
        :TrackID = mode(:TrackID),
        :Type = mode(:Type),
        :Gonad = mode(:Gonad),
        :Time = [:Time],
        :Duration = length(:Time),

        :CellPositionX = [:CellPositionX],
        :CellPositionY = [:CellPositionY],
        :CellPositionZ = [:CellPositionZ],

        :CellEllipsoidAxisCX = [:CellEllipsoidAxisCX],
        :CellEllipsoidAxisCY = [:CellEllipsoidAxisCY],
        :CellEllipsoidAxisCZ = [:CellEllipsoidAxisCZ],
        :MeanCellEllipsoidAxisCX = mean(:CellEllipsoidAxisCX),
        :MeanCellEllipsoidAxisCY = mean(:CellEllipsoidAxisCY),
        :MeanCellEllipsoidAxisCZ = mean(:CellEllipsoidAxisCZ),

        :CellEllipsoidAxisBX = [:CellEllipsoidAxisBX],
        :CellEllipsoidAxisBY = [:CellEllipsoidAxisBY],
        :CellEllipsoidAxisBZ = [:CellEllipsoidAxisBZ],
        :MeanCellEllipsoidAxisBX = mean(:CellEllipsoidAxisBX),
        :MeanCellEllipsoidAxisBY = mean(:CellEllipsoidAxisBY),
        :MeanCellEllipsoidAxisBZ = mean(:CellEllipsoidAxisBZ),

        :CellEllipsoidAxisAX = [:CellEllipsoidAxisAX],
        :CellEllipsoidAxisAY = [:CellEllipsoidAxisAY],
        :CellEllipsoidAxisAZ = [:CellEllipsoidAxisAZ],
        :MeanCellEllipsoidAxisAX = mean(:CellEllipsoidAxisAX),
        :MeanCellEllipsoidAxisAY = mean(:CellEllipsoidAxisAY),
        :MeanCellEllipsoidAxisAZ = mean(:CellEllipsoidAxisAZ),
        
        :CellEllipsoidAxisLengthC = [:CellEllipsoidAxisLengthC],
        :CellEllipsoidAxisLengthB = [:CellEllipsoidAxisLengthB],
        :CellEllipsoidAxisLengthA = [:CellEllipsoidAxisLengthA],



        :CellDistanceDTC = [:DistanceDTC],
        :MeanDistanceDTC = mean(:DistanceDTC),
        :CellVolume = [:CellVolume],
        :MeanCellVolume = mean(:CellVolume),
        :VarianceCellVolume = var(:CellVolume),
        :OrientationCtoDP = [:OrientationCtoDP],
        :MeanOrientationCtoDP = mean(:OrientationCtoDP),
        :VarianceOrientationCtoDP = var(:OrientationCtoDP),
        :ProlateIndex = [:Cell_Ellipticity_Prolate],
        :MeanProlateIndex = mean(:Cell_Ellipticity_Prolate))


    CellVolumeIntercept = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))
    CellVolumeSlope = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))
    CellVolumeR2 = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))
    CellVolumeStderrorIntercept = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))
    CellVolumeStderrorSloppe = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))
    CellSize = Vector{Union{String,Missing}}(undef, nrow(CellShape_df_final))
    CellOrientationIntercept = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))
    CellOrientationSlope = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))
    CellOrientationR2 = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))
    CellOrientationStderrorIntercept = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))
    CellOrientationStderrorSloppe = Vector{Union{Float64,Missing}}(undef, nrow(CellShape_df_final))

    for i in 1:nrow(CellShape_df_final)
        data = DataFrame(Time=CellShape_df_final[i, :Time], CellVolume=CellShape_df_final[i, :CellVolume], OrientationCtoDP=CellShape_df_final[i, :OrientationCtoDP])
        ols = lm(@formula(Time ~ CellVolume), data)
        ols2 = lm(@formula(Time ~ OrientationCtoDP), data)
        coeff_table = coeftable(ols)
        coefvalues = coef(ols)
        CellVolumeIntercept[i] = coefvalues[1]
        CellVolumeSlope[i] = coefvalues[2]
        CellVolumeR2[i] = r2(ols)
        CellVolumeStderrorIntercept[i] = stderror(ols)[1]
        CellVolumeStderrorSloppe[i] = stderror(ols)[2]
        if CellShape_df_final[i, :MeanCellVolume] <= 125
            CellSize[i] = "Small"
        else
            CellSize[i] = "Big"
        end
        coeff_table = coeftable(ols2)
        coefvalues = coef(ols2)
        CellOrientationIntercept[i] = coefvalues[1]
        CellOrientationSlope[i] = coefvalues[2]
        CellOrientationR2[i] = r2(ols2)
        CellOrientationStderrorIntercept[i] = stderror(ols2)[1]
        CellOrientationStderrorSloppe[i] = stderror(ols2)[2]
    end

    CellShape_df_final[!, :CellVolumeIntercept] = CellVolumeIntercept
    CellShape_df_final[!, :CellVolumeSlope] = CellVolumeSlope
    CellShape_df_final[!, :CellVolumeR2] = CellVolumeR2
    CellShape_df_final[!, :CellVolumeStderrorIntercept] = CellVolumeStderrorIntercept
    CellShape_df_final[!, :CellVolumeStderrorSloppe] = CellVolumeStderrorSloppe
    CellShape_df_final[!, :CellSize] = CellSize

    CellShape_df_final[!, :CellOrientationIntercept] = CellOrientationIntercept
    CellShape_df_final[!, :CellOrientationSlope] = CellOrientationSlope
    CellShape_df_final[!, :CellOrientationR2] = CellOrientationR2
    CellShape_df_final[!, :CellOrientationStderrorIntercept] = CellOrientationStderrorIntercept
    CellShape_df_final[!, :CellOrientationStderrorSloppe] = CellOrientationStderrorSloppe

    return CellShape_df_final
end


