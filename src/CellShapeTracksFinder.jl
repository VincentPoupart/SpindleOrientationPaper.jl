### Function to add cell shape data and mitotic steps to celloutput
function FindCellShapeTracksAndMitotiocSteps(Celloutput, CellShape, scoring_df)
    if hasproperty(Celloutput, :CellVolume) == false
        Celloutput[!, :TrackID] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :CellShapeTrackID] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellVolume] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :XPosition] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :YPosition] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ZPosition] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :NEBD] = Vector{Union{Missing,Int64}}(undef, nrow(Celloutput))
        Celloutput[!, :StartCongression] = Vector{Union{Missing,Int64}}(undef, nrow(Celloutput))
        Celloutput[!, :AnaphaseOnset] = Vector{Union{Missing,Int64}}(undef, nrow(Celloutput))
        Celloutput[!, :FrameFromAnaphaseOnset] = Vector{Vector{Union{Missing,Int64}}}(undef, nrow(Celloutput))
        Celloutput[!, :AngleLongAxisRachis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :AngleLongestAxisDP] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellVolume] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :CellVolumeSlope] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :VarianceCellVolume] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanAngleMeanCtoC] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :stdAngleMeanCtoC] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisCX] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisCY] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisCZ] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisBX] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisBY] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisBZ] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisAX] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisAY] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisAZ] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisLengthA] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisLengthB] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanCellEllipsoidAxisLengthC] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanDistanceDTC] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanOrientationCtoDP] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :VarianceOrientationCtoDP] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanProlateIndex] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :stdAngleMeanCtoC] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanAngleMeanCtoC] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanAngleSpindleRachisAnaphaseOnset] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanAngleSpindleLongAxisAnaphaseOnset] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanAngleLongAxisRachis] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :MeanAngleProjectLongestAxisDP] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisCX] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisCY] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisCZ] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisBX] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisBY] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisBZ] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisAX] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisAY] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisAZ] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProjectCellEllipsoidAxisCX] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProjectCellEllipsoidAxisCY] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProjectCellEllipsoidAxisCZ] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProjectCellEllipsoidAxisBX] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProjectCellEllipsoidAxisBY] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProjectCellEllipsoidAxisBZ] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProjectCellEllipsoidAxisAX] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProjectCellEllipsoidAxisAY] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProjectCellEllipsoidAxisAZ] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisLengthA] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisLengthB] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellEllipsoidAxisLengthC] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :CellDistanceDTC] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :OrientationCtoDP] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :ProlateIndex] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :AngleMeanCtoC] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        Celloutput[!, :AngleSpindleLongAxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
        for i in 1:nrow(Celloutput)
            Celloutput[i, :CellShapeTrackID] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellVolume] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :XPosition] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :YPosition] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ZPosition] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :FrameFromAnaphaseOnset] = Vector{Union{Missing,Int64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisCX] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisCY] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisCZ] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisBX] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisBY] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisBZ] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisAX] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisAY] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisAZ] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProjectCellEllipsoidAxisCX] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProjectCellEllipsoidAxisCY] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProjectCellEllipsoidAxisCZ] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProjectCellEllipsoidAxisBX] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProjectCellEllipsoidAxisBY] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProjectCellEllipsoidAxisBZ] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProjectCellEllipsoidAxisAX] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProjectCellEllipsoidAxisAY] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProjectCellEllipsoidAxisAZ] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisLengthA] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisLengthB] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellEllipsoidAxisLengthC] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :CellDistanceDTC] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :OrientationCtoDP] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :ProlateIndex] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :AngleMeanCtoC] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :AngleSpindleLongAxis] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :AngleLongAxisRachis] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            Celloutput[i, :AngleLongestAxisDP] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
            gonad = strip.(Celloutput[i, :gonad])
            cell = strip.(Celloutput[i, :cell])
            gonad_cell = string(gonad, " ", cell)
            IndexScoringDF = strip.(scoring_df[!, :Cell]) .== gonad_cell
            DPaxisVect = [missing, missing, missing]
            if haskey(DPaxis, gonad) == true
                DPaxisVect = [DPaxis[gonad][!, :X][1], DPaxis[gonad][!, :Y][1], DPaxis[gonad][!, :Z][1]]
            end
            if sum(IndexScoringDF) !== 0
                if isnan(scoring_df[IndexScoringDF, :NEBD][1]) == false
                    Celloutput[i, :NEBD] = scoring_df[IndexScoringDF, :NEBD][1]
                end
                if isnan(scoring_df[IndexScoringDF, :StartCongression][1]) == false
                    Celloutput[i, :StartCongression] = scoring_df[IndexScoringDF, :StartCongression][1]
                end
                if isnan(scoring_df[IndexScoringDF, :AnaphaseOnset][1]) == false
                    Celloutput[i, :AnaphaseOnset] = scoring_df[IndexScoringDF, :AnaphaseOnset][1]
                    Celloutput[i, :FrameFromAnaphaseOnset] = Celloutput[i, :Frame] .- Celloutput[i, :AnaphaseOnset]
                end
            end
            foo = CellShape[!, :Gonad] .== gonad
            for j in 1:Celloutput[i, :Duration]
                Celloutput[i, :XPosition][j] = mean([Celloutput[i, :CentaX][j], Celloutput[i, :CentbX][j]])
                Celloutput[i, :YPosition][j] = mean([Celloutput[i, :CentaY][j], Celloutput[i, :CentbY][j]])
                Celloutput[i, :ZPosition][j] = mean([Celloutput[i, :CentaZ][j], Celloutput[i, :CentbZ][j]])
                Frame = Celloutput[i, :Frame][j]
                poo = CellShape[foo, :]
                distances = Vector{Union{Missing,Float64}}(undef, nrow(poo))
                SpindleMidPoint = [Celloutput[i, :XPosition][j], Celloutput[i, :YPosition][j], Celloutput[i, :ZPosition][j]]
                for k in 1:nrow(poo)
                    FI = poo[k, :Time] .== Frame
                    if sum(FI) !== 0
                        CellPosition = [poo[k, :CellPositionX][FI][1], poo[k, :CellPositionY][FI][1], poo[k, :CellPositionZ][FI][1]]
                        distances[k] = Distance3D(SpindleMidPoint, CellPosition)
                    end
                end
                if isempty(skipmissing(distances)) == false
                    MinDistCellIndex = findmin(skipmissing(distances))
                    if MinDistCellIndex[1] <= 2.5
                        FI = poo[MinDistCellIndex[2], :Time] .== Frame
                        Celloutput[i, :CellShapeTrackID][j] = poo[MinDistCellIndex[2], :TrackID]
                        Celloutput[i, :CellVolume][j] = poo[MinDistCellIndex[2], :CellVolume][FI][1]
                        Celloutput[i, :CellEllipsoidAxisCX][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisCX][FI][1]
                        Celloutput[i, :CellEllipsoidAxisCY][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisCY][FI][1]
                        Celloutput[i, :CellEllipsoidAxisCZ][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisCZ][FI][1]
                        Celloutput[i, :CellEllipsoidAxisBX][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisBX][FI][1]
                        Celloutput[i, :CellEllipsoidAxisBY][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisBY][FI][1]
                        Celloutput[i, :CellEllipsoidAxisBZ][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisBZ][FI][1]
                        Celloutput[i, :CellEllipsoidAxisAX][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisAX][FI][1]
                        Celloutput[i, :CellEllipsoidAxisAY][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisAY][FI][1]
                        Celloutput[i, :CellEllipsoidAxisAZ][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisAZ][FI][1]
                        Celloutput[i, :CellEllipsoidAxisLengthA][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisLengthA][FI][1]
                        Celloutput[i, :CellEllipsoidAxisLengthB][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisLengthB][FI][1]
                        Celloutput[i, :CellEllipsoidAxisLengthC][j] = poo[MinDistCellIndex[2], :CellEllipsoidAxisLengthC][FI][1]
                        Celloutput[i, :CellDistanceDTC][j] = poo[MinDistCellIndex[2], :CellDistanceDTC][FI][1]
                        Celloutput[i, :OrientationCtoDP][j] = poo[MinDistCellIndex[2], :OrientationCtoDP][FI][1]
                        Celloutput[i, :ProlateIndex][j] = poo[MinDistCellIndex[2], :ProlateIndex][FI][1]
                    end
                end
            end
            if ismissing(Celloutput[i, :AnaphaseOnset]) == false
                IndexAnaphaseonset = findfirst(isequal(Celloutput[i, :AnaphaseOnset]), Celloutput[i, :Frame])
                if length(Celloutput[i, :Frame][1:IndexAnaphaseonset]) >= 5 && IndexAnaphaseonset + 1 <= length(Celloutput[i, :Frame])
                    data = DataFrame(Time=collect(0:1:4), CellVolume=Celloutput[i, :CellVolume][IndexAnaphaseonset-4:IndexAnaphaseonset])
                    println("i = ", i)
                    println(data)
                    if nrow(dropmissing(data)) == 5
                        ols = lm(@formula(Time ~ CellVolume), data)
                        coefvalues = coef(ols)
                        Celloutput[i, :MeanCellEllipsoidAxisCX] = mean(Celloutput[i, :CellEllipsoidAxisCX][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisCY] = mean(Celloutput[i, :CellEllipsoidAxisCY][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisCZ] = mean(Celloutput[i, :CellEllipsoidAxisCZ][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisBX] = mean(Celloutput[i, :CellEllipsoidAxisBX][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisBY] = mean(Celloutput[i, :CellEllipsoidAxisBY][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisBZ] = mean(Celloutput[i, :CellEllipsoidAxisBZ][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisAX] = mean(Celloutput[i, :CellEllipsoidAxisAX][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisAY] = mean(Celloutput[i, :CellEllipsoidAxisAY][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisAZ] = mean(Celloutput[i, :CellEllipsoidAxisAZ][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisLengthC] = mean(Celloutput[i, :CellEllipsoidAxisLengthC][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisLengthB] = mean(Celloutput[i, :CellEllipsoidAxisLengthB][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanCellEllipsoidAxisLengthA] = mean(Celloutput[i, :CellEllipsoidAxisLengthA][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        MeanC = [Celloutput[i, :MeanCellEllipsoidAxisCX], Celloutput[i, :MeanCellEllipsoidAxisCY], Celloutput[i, :MeanCellEllipsoidAxisCZ]]
                        for j in 1:Celloutput[i, :Duration]
                            LongAxis = [Celloutput[i, :CellEllipsoidAxisCX][j], Celloutput[i, :CellEllipsoidAxisCY][j], Celloutput[i, :CellEllipsoidAxisCZ][j]] .* Celloutput[i, :CellEllipsoidAxisLengthC][j]
                            MediumAxis = [Celloutput[i, :CellEllipsoidAxisBX][j], Celloutput[i, :CellEllipsoidAxisBY][j], Celloutput[i, :CellEllipsoidAxisBZ][j]] .* Celloutput[i, :CellEllipsoidAxisLengthB][j]
                            ShortAxis = [Celloutput[i, :CellEllipsoidAxisAX][j], Celloutput[i, :CellEllipsoidAxisAY][j], Celloutput[i, :CellEllipsoidAxisAZ][j]] .* Celloutput[i, :CellEllipsoidAxisLengthA][j]
                            Spindle = [Celloutput[i, :CentaX][j] - Celloutput[i, :CentbX][j], Celloutput[i, :CentaY][j] - Celloutput[i, :CentbY][j], Celloutput[i, :CentaZ][j] - Celloutput[i, :CentbZ][j]]
                            Rachis = [Celloutput[i, :XNormal][j], Celloutput[i, :YNormal][j], Celloutput[i, :ZNormal][j]]
                            if sum(ismissing.(LongAxis)) == 0 && sum(ismissing.(MeanC)) == 0
                                Celloutput[i, :AngleMeanCtoC][j] = CalculerAngle(LongAxis, MeanC)
                            end
                            if sum(ismissing.(LongAxis)) == 0 && sum(ismissing.(Spindle)) == 0
                                Celloutput[i, :AngleSpindleLongAxis][j] = CalculerAngle(LongAxis, Spindle)
                            end
                            if sum(ismissing.(LongAxis)) == 0 && sum(ismissing.(Rachis)) == 0
                                Celloutput[i, :AngleLongAxisRachis][j] = CalculerAngle(LongAxis, Rachis)
                                ProjectLongAxis = ProjectVonU(LongAxis, Rachis)
                                ProjectMediumAxis = ProjectVonU(MediumAxis, Rachis)
                                ProjectShortAxis = ProjectVonU(ShortAxis, Rachis)
                                Celloutput[i, :ProjectCellEllipsoidAxisCX][j] = ProjectLongAxis[1]
                                Celloutput[i, :ProjectCellEllipsoidAxisCY][j] = ProjectLongAxis[2]
                                Celloutput[i, :ProjectCellEllipsoidAxisCZ][j] = ProjectLongAxis[3]
                                Celloutput[i, :ProjectCellEllipsoidAxisBX][j] = ProjectMediumAxis[1]
                                Celloutput[i, :ProjectCellEllipsoidAxisBY][j] = ProjectMediumAxis[2]
                                Celloutput[i, :ProjectCellEllipsoidAxisBZ][j] = ProjectMediumAxis[3]
                                Celloutput[i, :ProjectCellEllipsoidAxisAX][j] = ProjectShortAxis[1]
                                Celloutput[i, :ProjectCellEllipsoidAxisAY][j] = ProjectShortAxis[2]
                                Celloutput[i, :ProjectCellEllipsoidAxisAZ][j] = ProjectShortAxis[3]
                                ProjectDPaxisVect = ProjectVonU(DPaxisVect, Rachis)
                                NormsProjectAxes = Vector{Union{Missing,Float64}}(undef, 3)
                                NormsProjectAxes[1] = norm(ProjectLongAxis)
                                NormsProjectAxes[2] = norm(ProjectMediumAxis)
                                NormsProjectAxes[3] = norm(ProjectShortAxis)
                                IndexLongerProjectAxis = findmax(skipmissing(NormsProjectAxes))
                                LongestProjectAxis = []
                                if IndexLongerProjectAxis[2] == 1
                                    LongestProjectAxis = ProjectLongAxis
                                elseif IndexLongerProjectAxis[2] == 2
                                    LongestProjectAxis = ProjectMediumAxis
                                else
                                    LongestProjectAxis = ProjectShortAxis
                                end
                                Celloutput[i, :AngleLongestAxisDP][j] = CalculerAngle(LongestProjectAxis, ProjectDPaxisVect)
                            end
                        end
                        Celloutput[i, :TrackID] = mode(skipmissing(Celloutput[i, :CellShapeTrackID]))
                        Celloutput[i, :MeanCellVolume] = mean(Celloutput[i, :CellVolume][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :CellVolumeSlope] = coefvalues[2]
                        Celloutput[i, :VarianceCellVolume] = var(Celloutput[i, :CellVolume][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanAngleMeanCtoC] = mean(Celloutput[i, :AngleMeanCtoC][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :stdAngleMeanCtoC] = std(Celloutput[i, :AngleMeanCtoC][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanAngleSpindleRachisAnaphaseOnset] = mean(Celloutput[i, :AngleSpindleRachis][IndexAnaphaseonset-1:IndexAnaphaseonset+1])
                        Celloutput[i, :MeanDistanceDTC] = mean(Celloutput[i, :CellDistanceDTC][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanOrientationCtoDP] = mean(Celloutput[i, :OrientationCtoDP][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :VarianceOrientationCtoDP] = var(Celloutput[i, :OrientationCtoDP][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanProlateIndex] = mean(Celloutput[i, :ProlateIndex][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :stdAngleMeanCtoC] = std(Celloutput[i, :AngleMeanCtoC][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanAngleMeanCtoC] = mean(Celloutput[i, :AngleMeanCtoC][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanAngleSpindleLongAxisAnaphaseOnset] = mean(Celloutput[i, :AngleSpindleLongAxis][IndexAnaphaseonset-1:IndexAnaphaseonset+1])
                        Celloutput[i, :MeanAngleLongAxisRachis] = mean(Celloutput[i, :AngleLongAxisRachis][IndexAnaphaseonset-4:IndexAnaphaseonset])
                        Celloutput[i, :MeanAngleProjectLongestAxisDP] = mean(Celloutput[i, :AngleLongestAxisDP][IndexAnaphaseonset-4:IndexAnaphaseonset])
                    end
                end
            end

        end
    end
    return Celloutput
end

