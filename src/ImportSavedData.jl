#### Set the path to saved data
path4 = "M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/"

### Import saved scoring data (NEBD, Start Congression, Anaphase Onset)
scoring_df = CSV.read(joinpath(path2, "scoring_df.csv"), DataFrame)

### Import saved Interphase data
Interphase_df = DataFrame(CSV.File(joinpath(path4, "InterphaseDF.csv"), delim=","))
Interphase_df = ParseDataframe(Interphase_df)
Interphase_df = dropmissing(Interphase_df, :MeanOrientationCtoDP)

### Import saved Mitosis data
Mitosis_df = DataFrame(CSV.File(joinpath(path4, "MitosisDF.csv"), delim=","))
Mitosis_df = ParseDataframe(Mitosis_df)
Mitosis_df = dropmissing(Mitosis_df, [:MeanOrientationCtoDP, :MeanAngleLongAxisRachis])


function Sure10TimePointsInterphase(df)
    namesdf = names(df)
    for i in 1:nrow(df)
        for name in namesdf
            if name == "MeanCellEllipsoidAxisCX"
                df[i, name] = mean(df[i, :CellEllipsoidAxisCX][1:10])
            elseif  name == "MeanCellEllipsoidAxisCY"
                df[i, name] = mean(df[i, :CellEllipsoidAxisCY][1:10])
            elseif  name == "MeanCellEllipsoidAxisCZ"
                df[i, name] = mean(df[i, :CellEllipsoidAxisCZ][1:10])
            elseif  name == "MeanCellEllipsoidAxisBX"
                df[i, name] = mean(df[i, :CellEllipsoidAxisBX][1:10])
            elseif  name == "MeanCellEllipsoidAxisBY"
                df[i, name] = mean(df[i, :CellEllipsoidAxisBY][1:10])
            elseif  name == "MeanCellEllipsoidAxisBZ"
                df[i, name] = mean(df[i, :CellEllipsoidAxisBZ][1:10])
            elseif  name == "MeanCellEllipsoidAxisAX"
                df[i, name] = mean(df[i, :CellEllipsoidAxisAX][1:10])
            elseif  name == "MeanCellEllipsoidAxisAY"
                df[i, name] = mean(df[i, :CellEllipsoidAxisAY][1:10])
            elseif  name == "MeanCellEllipsoidAxisAZ"
                df[i, name] = mean(df[i, :CellEllipsoidAxisAZ][1:10])
            elseif  name == "MeanDistanceDTC"
                df[i, name] = mean(df[i, :CellDistanceDTC][1:10])  
            elseif  name == "MeanCellVolume"
                df[i, name] = mean(df[i, :CellVolume][1:10])
            elseif  name == "MeanOrientationCtoDP"
                df[i, name] = mean(df[i, :OrientationCtoDP][1:10])
            elseif  name == "MeanProlateIndex"
                df[i, name] = mean(df[i, :ProlateIndex][1:10])
            elseif  name == "MeanAnglesLongAxisRachis"
                df[i, name] = mean(df[i, :AngleLongAxisRachis][1:10])
            elseif  name == "MeanAnglesProjectLongestAxisDP"
                df[i, name] = mean(df[i, :AnglesProjectLongestAxisDP][1:10])
            
            end
        end
    end
    return df
end


Interphase_df = Sure10TimePointsInterphase(Interphase_df)