####### Function to add surface normal vector

function NormalVector2MitosisDict!(Mitosis_dict, wrlfilelist; algorythm="GrowingSphere")
    Mitosis_dict = sort(Mitosis_dict)
    k = sort(collect(keys(Mitosis_dict)))
    wrlbasenamelist = RemoveExt.(basename.(wrlfilelist), ".")
    GonadWRL = nothing
    IndexedFaceSets = nothing
    for i in eachindex(k)
        foo = []
        Gonad_df = RemoveExt(k[i], " Cell")
        foo = Gonad_df .== wrlbasenamelist
        if hasproperty(Mitosis_dict[k[i]], :AngleSpindleRachis) == false && sum(foo) > 0
            if isnothing(IndexedFaceSets) || i == 1 || Gonad_df !== RemoveExt(k[i-1], " Cell")
                IndexedFaceSets = nothing
                GonadWRL = nothing
                wrl = unique(wrlfilelist[foo])
                GonadWRL = RemoveExt.(basename.(wrl[1]), ".")
                IndexedFaceSets = ExtractDataVRML(wrl[1])
                println("gonad ", Gonad_df, " is done")
            end
            if Gonad_df == GonadWRL
                if algorythm == "GrowingSphere"
                    Mitosis_dict[k[i]] = AddAnglesGrowingSphereMitosis(Mitosis_dict[k[i]], IndexedFaceSets)
                elseif algorythm == "First50um"
                    Mitosis_dict[k[i]] = AddAnglesFirst50umMitosis(Mitosis_dict[k[i]], IndexedFaceSets)
                end
            end
        end
        println("cell $i is done")
    end
    return Mitosis_dict
end


function AddAnglesFirst50umMitosis(fazer::DataFrame, IndexedFaceSets)
    NormalsX = Vector{Union{Float64,Missing}}(undef, length(fazer[!, :FRAME]))
    NormalsY = Vector{Union{Float64,Missing}}(undef, length(fazer[!, :FRAME]))
    NormalsZ = Vector{Union{Float64,Missing}}(undef, length(fazer[!, :FRAME]))
    Angles = Vector{Union{Float64,Missing}}(undef, length(fazer[!, :FRAME]))
    for j in 1:nrow(fazer)
        Frame = fazer[j, :FRAME]
        if Frame <= length(IndexedFaceSets)
            spindle_mid_point = fazer[j, :Spindle_Midpoint]
            cent = IndexedFaceSets[Frame]["centroids"]
            d = Vector{Float64}(undef, 0)
            for m in eachindex(cent)
                push!(d, Distance3D(cent[m], spindle_mid_point))
            end
            Ioo = sortperm(d)
            Sorted_Normals = IndexedFaceSets[Frame]["normals"][Ioo]
            Sorted_Areas = IndexedFaceSets[Frame]["areas"][Ioo]
            CumSumArea = cumsum(Sorted_Areas)
            IndexArea = findfirst(x -> x >= 50.0, CumSumArea)
            if IndexArea !== nothing && d[Ioo][IndexArea] <= 8
                spindle = [fazer[j, :POSITION_X] - fazer[j, :POSITION_X_1] fazer[j, :POSITION_Y] - fazer[j, :POSITION_Y_1] fazer[j, :POSITION_Z] - fazer[j, :POSITION_Z_1]]
                normal = sum(Sorted_Normals[1:IndexArea])
                NormalsX[j] = normal[1]
                NormalsY[j] = normal[2]
                NormalsZ[j] = normal[3]
                Angles[j] = CalculerAngle(normal, spindle)
            end
        end
    end
    fazer[!, :XNormal] = NormalsX
    fazer[!, :YNormal] = NormalsY
    fazer[!, :ZNormal] = NormalsZ
    fazer[!, :AngleSpindleRachis] = Angles
    return fazer
end

function AddAnglesGrowingSphereMitosis(fazer::DataFrame, IndexedFaceSets)
    fazer[!, :XNormal] = Vector{Union{Float64,Missing}}(undef, length(fazer[!, :FRAME]))
    fazer[!, :YNormal] = Vector{Union{Float64,Missing}}(undef, length(fazer[!, :FRAME]))
    fazer[!, :ZNormal] = Vector{Union{Float64,Missing}}(undef, length(fazer[!, :FRAME]))
    Angles = Vector{Union{Float64,Missing}}(undef, length(fazer[!, :FRAME]))
    for j in 1:nrow(fazer)
        Frame = fazer[j, :FRAME]
        if Frame <= length(IndexedFaceSets)
            spindle_mid_point = fazer[j, :Spindle_Midpoint]
            cent = IndexedFaceSets[Frame]["centroids"]
            d = Vector{Union{Missing,Float64}}(undef, length(cent))
            for m in eachindex(cent)
                d[m] = Distance3D(cent[m], spindle_mid_point)
            end
            Ioo = sortperm(d)
            Sorted_Normals = IndexedFaceSets[Frame]["normals"][Ioo]
            Sorted_Areas = IndexedFaceSets[Frame]["areas"][Ioo]
            CumSumArea = cumsum(Sorted_Areas)
            distance = 1
            area = 0
            IndexDistance = 0
            if last(CumSumArea) > 100
                while area < 50
                    IndexDistance = findfirst(x -> x >= distance, d[Ioo])
                    area = CumSumArea[IndexDistance]
                    distance += 0.1
                end
                println("area = ", area)
                if IndexDistance !== nothing && d[Ioo][IndexDistance] <= 8
                    spindle = [fazer[j, :POSITION_X] - fazer[j, :POSITION_X_1] fazer[j, :POSITION_Y] - fazer[j, :POSITION_Y_1] fazer[j, :POSITION_Z] - fazer[j, :POSITION_Z_1]]
                    normal = sum(Sorted_Normals[1:IndexDistance])
                    fazer[j, :XNormal] = normal[1]
                    fazer[j, :YNormal] = normal[2]
                    fazer[j, :ZNormal] = normal[3]
                    Angles[j] = CalculerAngle(normal, spindle)
                end
            end
        end
    end
    fazer[!, :AngleSpindleRachis] = Angles
    return fazer
end

function NormalVector2InterphaseDF!(Interphase_df, wrlfilelist; algorythm="GrowingSphere")
    wrlbasenamelist = RemoveExt.(basename.(wrlfilelist), ".")
    GonadWRL = nothing
    IndexedFaceSets = nothing
    if hasproperty(Interphase_df, :AngleLongAxisRachis) == false
        Interphase_df[!, :AngleLongAxisRachis] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Interphase_df))
        Interphase_df[!, :XNormal] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Interphase_df))
        Interphase_df[!, :YNormal] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Interphase_df))
        Interphase_df[!, :ZNormal] = Vector{Vector{Union{Missing,Float64}}}(undef, nrow(Interphase_df))
        for i in 1:nrow(Interphase_df)
            Interphase_df[i, :AngleLongAxisRachis] = Vector{Union{Missing,Float64}}(undef, Interphase_df[i, :Duration])
            Interphase_df[i, :XNormal] = Vector{Union{Missing,Float64}}(undef, Interphase_df[i, :Duration])
            Interphase_df[i, :YNormal] = Vector{Union{Missing,Float64}}(undef, Interphase_df[i, :Duration])
            Interphase_df[i, :ZNormal] = Vector{Union{Missing,Float64}}(undef, Interphase_df[i, :Duration])
        end
    end
    for i in 1:nrow(Interphase_df)
        foo = []
        Gonad_df = Interphase_df[i, :Gonad]
        foo = Gonad_df .== wrlbasenamelist
        if sum(ismissing.(Interphase_df[i, :AngleLongAxisRachis]) .== false) == 0 && sum(foo) > 0
            if isnothing(IndexedFaceSets) || i == 1 || Gonad_df !== Interphase_df[i-1, :Gonad]
                IndexedFaceSets = nothing
                GonadWRL = nothing
                wrl = unique(wrlfilelist[foo])
                GonadWRL = RemoveExt.(basename.(wrl[1]), ".")
                IndexedFaceSets = ExtractDataVRML(wrl[1])
                println("gonad ", Gonad_df, " is done")
            end
            if Gonad_df == GonadWRL
                if algorythm == "GrowingSphere"
                    Interphase_df[i, :] = AddAnglesGrowingSphereMitosis!(Interphase_df[i, :], IndexedFaceSets)
                elseif algorythm == "First50um"
                    Interphase_df[i, :] = AddAnglesFirst50umInterphase!(Interphase_df[i, :], IndexedFaceSets)
                end
            end
        end
        println("cell $i is done")
    end
    return Interphase_df
end

function AddAnglesFirst50umInterphase!(fazer::DataFrameRow{DataFrame,DataFrames.Index}, IndexedFaceSets)
    NormalsX = Vector{Union{Missing,Float64}}(undef, fazer[:Duration])
    NormalsY = Vector{Union{Missing,Float64}}(undef, fazer[:Duration])
    NormalsZ = Vector{Union{Missing,Float64}}(undef, fazer[:Duration])
    Angles = Vector{Union{Missing,Float64}}(undef, fazer[:Duration])
    for j in 1:fazer[:Duration]
        Frame = fazer[:Time][j]
        if Frame <= length(IndexedFaceSets)
            spindle_mid_point = [fazer[:CellPositionX][j] fazer[:CellPositionY][j] fazer[:CellPositionZ][j]]
            cent = IndexedFaceSets[Frame]["centroids"]
            d = Vector{Union{Missing,Float64}}(undef, length(cent))
            for m in eachindex(cent)
                d[m] = Distance3D(cent[m], spindle_mid_point)
            end
            Ioo = sortperm(d)
            Sorted_Normals = IndexedFaceSets[Frame]["normals"][Ioo]
            Sorted_Areas = IndexedFaceSets[Frame]["areas"][Ioo]
            CumSumArea = cumsum(Sorted_Areas)
            if length(CumSumArea) > 100
                IndexArea = findfirst(x -> x >= 50.0, CumSumArea)
                if IndexArea !== nothing && d[Ioo][IndexArea] <= 8 && last(CumSumArea) > 100
                    LongAxis = [fazer[:CellEllipsoidAxisCX][j] fazer[:CellEllipsoidAxisCY][j] fazer[:CellEllipsoidAxisCZ][j]]
                    normal = sum(Sorted_Normals[1:IndexArea])
                    NormalsX[j] = normal[1]
                    NormalsY[j] = normal[2]
                    NormalsZ[j] = normal[3]
                    Angles[j] = CalculerAngle(normal, LongAxis)
                end
            end
        end
    end

    fazer[:XNormal] = NormalsX
    fazer[:YNormal] = NormalsY
    fazer[:ZNormal] = NormalsZ
    fazer[:AngleLongAxisRachis] = Angles
    return fazer
end

function AddAnglesGrowingSphereInterphase!(fazer::DataFrameRow{DataFrame,DataFrames.Index}, IndexedFaceSets)
    NormalsX = Vector{Union{Missing,Float64}}(undef, fazer[:Duration])
    NormalsY = Vector{Union{Missing,Float64}}(undef, fazer[:Duration])
    NormalsZ = Vector{Union{Missing,Float64}}(undef, fazer[:Duration])
    Angles = Vector{Union{Missing,Float64}}(undef, fazer[:Duration])
    for j in 1:fazer[:Duration]
        Frame = fazer[:Time][j]
        if Frame <= length(IndexedFaceSets)
            spindle_mid_point = [fazer[:CellPositionX][j] fazer[:CellPositionY][j] fazer[:CellPositionZ][j]]
            cent = IndexedFaceSets[Frame]["centroids"]
            d = Vector{Union{Missing,Float64}}(undef, length(cent))
            for m in eachindex(cent)
                d[m] = Distance3D(cent[m], spindle_mid_point)
            end
            Ioo = sortperm(d)
            Sorted_Normals = IndexedFaceSets[Frame]["normals"][Ioo]
            Sorted_Areas = IndexedFaceSets[Frame]["areas"][Ioo]
            CumSumArea = cumsum(Sorted_Areas)
            CumSumArea = cumsum(Sorted_Areas)
            distance = 1
            area = 0
            IndexDistance = 0
            if last(CumSumArea) > 100
                while area < 50
                    IndexDistance = findfirst(x -> x >= distance, d[Ioo])
                    area = CumSumArea[IndexDistance]
                    distance += 0.1
                end
                println("area = ", area)
                if IndexDistance !== nothing && d[Ioo][IndexDistance] <= 8
                    LongAxis = [fazer[:CellEllipsoidAxisCX][j] fazer[:CellEllipsoidAxisCY][j] fazer[:CellEllipsoidAxisCZ][j]]
                    normal = sum(Sorted_Normals[1:IndexDistance])
                    NormalsX[j] = normal[1]
                    NormalsY[j] = normal[2]
                    NormalsZ[j] = normal[3]
                    Angles[j] = CalculerAngle(normal, LongAxis)
                end
            end
        end
    end
    fazer[:XNormal] = NormalsX
    fazer[:YNormal] = NormalsY
    fazer[:ZNormal] = NormalsZ
    fazer[:AngleLongAxisRachis] = Angles
    return fazer
end


Mitosis_dict = NormalVector2MitosisDF!(Mitosis_dict, wrlfiles; algorythm="GrowingSphere")
