function AddProjectAxis(CellShape, DPaxis)
    CellShape[!, :MeanAngleProjectLongestAxisDP] = Vector{Union{Missing,Float64}}(undef, nrow(CellShape))
    CellShape[!, :MeanAnglesLongAxisRachis] = Vector{Union{Missing,Float64}}(undef, nrow(CellShape))

    CellShape[!, :AngleProjectLongestAxisDP] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(CellShape))

    for i in 1:nrow(CellShape)

        CellShape[i, :AngleProjectLongestAxisDP] = Vector{Union{Missing,Float64}}(undef, CellShape[i, :Duration])
        gonad = strip.(CellShape[i, :Gonad])
        DPaxisVect = [missing, missing, missing]
        if haskey(DPaxis, gonad) == true
            DPaxisVect = [DPaxis[gonad][!, :X][1], DPaxis[gonad][!, :Y][1], DPaxis[gonad][!, :Z][1]]
            println("i = ", i)
            for j in 1:CellShape[i, :Duration]
                LongAxis = [CellShape[i, :CellEllipsoidAxisCX][j], CellShape[i, :CellEllipsoidAxisCY][j], CellShape[i, :CellEllipsoidAxisCZ][j]] .* CellShape[i, :CellEllipsoidAxisLengthC][j]
                MediumAxis = [CellShape[i, :CellEllipsoidAxisBX][j], CellShape[i, :CellEllipsoidAxisBY][j], CellShape[i, :CellEllipsoidAxisBZ][j]] .* CellShape[i, :CellEllipsoidAxisLengthB][j]
                ShortAxis = [CellShape[i, :CellEllipsoidAxisAX][j], CellShape[i, :CellEllipsoidAxisAY][j], CellShape[i, :CellEllipsoidAxisAZ][j]] .* CellShape[i, :CellEllipsoidAxisLengthA][j]
                Rachis = [CellShape[i, :XNormal][j], CellShape[i, :YNormal][j], CellShape[i, :ZNormal][j]]
                if sum(ismissing.(Rachis)) == 0 && sum(ismissing.(LongAxis)) == 0
                    ProjectLongAxis = ProjectVonU(LongAxis, Rachis)
                    ProjectMediumAxis = ProjectVonU(MediumAxis, Rachis)
                    ProjectShortAxis = ProjectVonU(ShortAxis, Rachis)
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
                    CellShape[i, :AngleProjectLongestAxisDP][j] = CalculerAngle(LongestProjectAxis, ProjectDPaxisVect)
                end
            end
        end
        CellShape[i, :MeanAngleProjectLongestAxisDP] = mean(CellShape[i, :AngleProjectLongestAxisDP][1:10])
        CellShape[i, :MeanAnglesLongAxisRachis] = mean(CellShape[i, :AngleLongAxisRachis][1:10])
    end
    return CellShape
end

function AddProjectSpindleAngleCelloutput(Celloutput, DPaxis)
    Celloutput[!, :MeanAngleProjectSpindleProjectLongAxis] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
    Celloutput[!, :AngleProjectSpindleProjectLongAxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))

    Celloutput[!, :MeanAngleProjectSpindleProjectShortAxis] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
    Celloutput[!, :AngleProjectSpindleProjectShortAxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))

    Celloutput[!, :MeanAngleProjectSpindleProjectDPaxis] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
    Celloutput[!, :AngleProjectSpindleProjectDPaxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))

    Celloutput[!, :MeanAngleProjectDPaxisProjectLongAxis] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
    Celloutput[!, :AngleProjectDPaxisProjectLongAxis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))
    
    Celloutput[!, :MeanLongestProjectedAxisNorm] = Vector{Union{Missing,Float64}}(undef, nrow(Celloutput))
    Celloutput[!, :LongestProjectedAxisNorm] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Celloutput))


    
    for i in 1:nrow(Celloutput)
        
        if ismissing(Celloutput[i, :AnaphaseOnset]) == false
            
            IndexAnaphaseonset = findfirst(isequal(Celloutput[i, :AnaphaseOnset]), Celloutput[i, :Frame])
            if length(Celloutput[i, :Frame][1:IndexAnaphaseonset]) >= 5 && IndexAnaphaseonset + 1 <= length(Celloutput[i, :Frame])
                Celloutput[i, :AngleProjectSpindleProjectLongAxis] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
                Celloutput[i, :AngleProjectSpindleProjectShortAxis] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
                Celloutput[i, :AngleProjectSpindleProjectDPaxis] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
                Celloutput[i, :AngleProjectDPaxisProjectLongAxis] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])
                Celloutput[i, :LongestProjectedAxisNorm] = Vector{Union{Missing,Float64}}(undef, Celloutput[i, :Duration])

                gonad = strip.(Celloutput[i, :gonad])
                DPaxisVect = [missing, missing, missing]
            
                if haskey(DPaxis, gonad) == true
                    DPaxisVect = [DPaxis[gonad][!, :X][1], DPaxis[gonad][!, :Y][1], DPaxis[gonad][!, :Z][1]]
                    for j in 1:Celloutput[i, :Duration]
                        LongAxis = [Celloutput[i, :CellEllipsoidAxisCX][j], Celloutput[i, :CellEllipsoidAxisCY][j], Celloutput[i, :CellEllipsoidAxisCZ][j]] .* Celloutput[i, :CellEllipsoidAxisLengthC][j]
                        MediumAxis = [Celloutput[i, :CellEllipsoidAxisBX][j], Celloutput[i, :CellEllipsoidAxisBY][j], Celloutput[i, :CellEllipsoidAxisBZ][j]] .* Celloutput[i, :CellEllipsoidAxisLengthB][j]
                        ShortAxis = [Celloutput[i, :CellEllipsoidAxisAX][j], Celloutput[i, :CellEllipsoidAxisAY][j], Celloutput[i, :CellEllipsoidAxisAZ][j]] .* Celloutput[i, :CellEllipsoidAxisLengthA][j]
                        Rachis = [Celloutput[i, :XNormal][j], Celloutput[i, :YNormal][j], Celloutput[i, :ZNormal][j]]
                        Spindle = [Celloutput[i, :CentaX][j] - Celloutput[i, :CentbX][j], Celloutput[i, :CentaY][j] - Celloutput[i, :CentbY][j], Celloutput[i, :CentaZ][j] - Celloutput[i, :CentbZ][j]]
                        if sum(ismissing.(Rachis)) == 0 && sum(ismissing.(LongAxis)) == 0
                            ProjectLongAxis = ProjectVonU(LongAxis, Rachis)
                            ProjectMediumAxis = ProjectVonU(MediumAxis, Rachis)
                            ProjectShortAxis = ProjectVonU(ShortAxis, Rachis)
                            ProjectSpindle = ProjectVonU(Spindle, Rachis)
                            ProjectDPaxisVect = ProjectVonU(DPaxisVect, Rachis)

                            NormsProjectAxes = Vector{Union{Missing,Float64}}(undef, 3)
                            NormsProjectAxes[1] = norm(ProjectLongAxis)
                            NormsProjectAxes[2] = norm(ProjectMediumAxis)
                            NormsProjectAxes[3] = norm(ProjectShortAxis)

                            IndexLongerProjectAxis = findmax(skipmissing(NormsProjectAxes))
                            IndexShorterProjectAxis = findmin(skipmissing(NormsProjectAxes))

                            LongestProjectAxis = []
                            ShortestProjectAxis = []
                            if IndexLongerProjectAxis[2] == 1
                                LongestProjectAxis = ProjectLongAxis
                            elseif IndexLongerProjectAxis[2] == 2
                                LongestProjectAxis = ProjectMediumAxis
                            else
                                LongestProjectAxis = ProjectShortAxis
                            end

                            if IndexShorterProjectAxis[2] == 1
                                ShortestProjectAxis = ProjectLongAxis
                            elseif IndexShorterProjectAxis[2] == 2
                                ShortestProjectAxis = ProjectMediumAxis
                            else
                                ShortestProjectAxis = ProjectShortAxis
                            end

                            Celloutput[i, :AngleProjectSpindleProjectLongAxis][j] = CalculerAngle(LongestProjectAxis, ProjectSpindle)
                            Celloutput[i, :AngleProjectSpindleProjectShortAxis][j] = CalculerAngle(ShortestProjectAxis, ProjectSpindle)
                            Celloutput[i, :AngleProjectSpindleProjectDPaxis][j] = CalculerAngle(ProjectDPaxisVect, ProjectSpindle)
                            Celloutput[i, :AngleProjectDPaxisProjectLongAxis][j] = CalculerAngle(ProjectDPaxisVect, LongestProjectAxis)
                            Celloutput[i, :LongestProjectedAxisNorm][j] = norm(LongestProjectAxis)
                            println(LongestProjectAxis)
                        end
                    end
                end
            end
        end
        Celloutput[i, :MeanAngleProjectSpindleProjectLongAxis] = mean(Celloutput[i, :AngleProjectSpindleProjectLongAxis][IndexAnaphaseonset-4:IndexAnaphaseonset])
        Celloutput[i, :MeanAngleProjectSpindleProjectShortAxis] = mean(Celloutput[i, :AngleProjectSpindleProjectShortAxis][IndexAnaphaseonset-4:IndexAnaphaseonset])
        Celloutput[i, :MeanAngleProjectSpindleProjectDPaxis] = mean(Celloutput[i, :AngleProjectSpindleProjectDPaxis][IndexAnaphaseonset-4:IndexAnaphaseonset])
        Celloutput[i, :MeanAngleProjectDPaxisProjectLongAxis] = mean(Celloutput[i, :AngleProjectDPaxisProjectLongAxis][IndexAnaphaseonset-4:IndexAnaphaseonset])
        Celloutput[i, :MeanLongestProjectedAxisNorm] = mean(Celloutput[i, :LongestProjectedAxisNorm][IndexAnaphaseonset-4:IndexAnaphaseonset])
    
    end
    return Celloutput
end

