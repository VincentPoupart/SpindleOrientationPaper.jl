

CairoMakie.activate!(type = "png")
colors = Makie.wong_colors()

filtered_df_big = filter(row -> row[:CellSize] == "Big", Interphase_df)
filtered_df_small = filter(row -> row[:CellSize] == "Small", Interphase_df)
Interphase_big_prolate = filter(row -> row[:MeanProlateIndex] >= 0.4, Interphase_df)
Interphase_small_prolate = filter(row -> row[:MeanProlateIndex] < 0.4, Interphase_df)


#### Theretical random polar histogram
n = 10000000
ys = rand(n)
xs = acos.(1 .- ys) .* (180 / π)
f1 = CairoMakie.Figure(size = (600,600))
PolarHist(xs, f1, "Theoretical random polar histogram", 2, 0.75, 1)
f1

f14 = CairoMakie.Figure(size = (1350, 450))
ax14 = CairoMakie.Axis(f14[1,1], 
    title = " Cells long axis orientation \n at interphase, n = $(nrow(Interphase_df))", 
    xlabel = "Cells long axis orientation (degrees)", 
    ylabel = "fraction", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)
DP = CairoMakie.stairs!(ax14, sort(Interphase_df[!, :MeanOrientationCtoDP]), ecdf(sort(Interphase_df[!, :MeanOrientationCtoDP]))(sort(Interphase_df[!, :MeanOrientationCtoDP])), color = Makie.wong_colors()[1])
rachis = CairoMakie.stairs!(ax14, sort(Interphase_df[!, :MeanAnglesLongAxisRachis]), ecdf(sort(Interphase_df[!, :MeanAnglesLongAxisRachis]))(sort(Interphase_df[!, :MeanAnglesLongAxisRachis])), color = Makie.wong_colors()[3])
random = CairoMakie.lines!(ax14, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1, color = Makie.wong_colors()[7])
PolarHist(Interphase_df[!, :MeanOrientationCtoDP], f14, "D/P axis", 1, 0.75, 3)
PolarHist(Interphase_df[!, :MeanAnglesLongAxisRachis], f14, "Rachis surface normal", 3, 0.75, 4)
CairoMakie.Legend(f14[1,2], [DP, rachis, random], ["D/P axis", "rachis surface normal", "theoretical random"])
CairoMakie.display(f14)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/CellsLongAxisOrientationInterphase1.pdf", f14, px_per_unit = 2) 


f16 = Figure(size = (700, 450))
ax16 = Axis(f16[1,1], 
    title = "Cell volume (µm³) according to frame from anaphse onset, n = $(nrow(Mitosis_df))", 
    xlabel = "Frame from anaphase onset", 
    ylabel = "Cell volume (µm³)",
    limits = (-5, 10, 50, 300))
for i in 1:nrow(Mitosis_df)
    println(i)
    lines!(ax16, Mitosis_df[i, :FrameFromAnaphaseOnset], Mitosis_df[i, :CellVolume])
end
display(f16)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/CellsVolumeFrameFromAna.pdf", f16, px_per_unit = 2)

f17 = CairoMakie.Figure(size = (2300, 450))
ax17 = CairoMakie.Axis(f17[1,1], 
    title = " Cells mitotic spindle orientation \n at anaphase onset, n = $(nrow(Mitosis_df))", 
    xlabel = "Mitotic spindle orientation (degrees)", 
    ylabel = "fraction", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)
ShortAxis = CairoMakie.stairs!(ax17, sort(Mitosis_df[!, :MeanAngleSpindleShortAxis]), ecdf(sort(Mitosis_df[!, :MeanAngleSpindleShortAxis]))(sort(Mitosis_df[!, :MeanAngleSpindleShortAxis])), color = (Makie.wong_colors()[5]))
LongAxis = CairoMakie.stairs!(ax17, sort(Mitosis_df[!, :MeanAngleSpindleLongAxisAnaphaseOnset]), ecdf(sort(Mitosis_df[!, :MeanAngleSpindleLongAxisAnaphaseOnset]))(sort(Mitosis_df[!, :MeanAngleSpindleLongAxisAnaphaseOnset])), color = (Makie.wong_colors()[1]))
rachis = CairoMakie.stairs!(ax17, sort(Mitosis_df[!, :MeanAngleSpindleRachisAnaphaseOnset]), ecdf(sort(Mitosis_df[!, :MeanAngleSpindleRachisAnaphaseOnset]))(sort(Mitosis_df[!, :MeanAngleSpindleRachisAnaphaseOnset])), color = Makie.wong_colors()[3])
DP =  CairoMakie.stairs!(ax17, sort(Mitosis_df[!, :MeanAngleSpindleDPaxis]), ecdf(sort(Mitosis_df[!, :MeanAngleSpindleDPaxis]))(sort(Mitosis_df[!, :MeanAngleSpindleDPaxis])), color = Makie.wong_colors()[4])
random = CairoMakie.lines!(ax17, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1, color = Makie.wong_colors()[2])
CairoMakie.Legend(f17[1,2], [ShortAxis, LongAxis, rachis, DP, random], ["cell short axis","cell long axis", "rachis surface normal ","D/P axis", "theoretical random"])
PolarHist(Mitosis_df[!, :MeanAngleSpindleShortAxis], f17, "Cell short axis", 5, 0.75, 3)
PolarHist(Mitosis_df[!, :MeanAngleSpindleLongAxisAnaphaseOnset], f17, "Cell long axis", 1, 0.75, 4)
PolarHist(Mitosis_df[!, :MeanAngleSpindleRachisAnaphaseOnset], f17, "Rachis surface normal", 3, 0.75, 5)
PolarHist(Mitosis_df[!, :MeanAngleSpindleDPaxis], f17, "D/P axis", 4, 0.75, 6)
CairoMakie.display(f17)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/CellsSpindleOrientationMitosis2.pdf", f17, px_per_unit = 2) 


f18 = CairoMakie.Figure(size = (2250, 450))
ax18 = CairoMakie.Axis(f18[1,1], 
    title = " Cells long axis orientation to D/P axis and rachis normal \n at interphase, n (small) = $(nrow(filtered_df_small)), n (big) = $(nrow(filtered_df_big))", 
    xlabel = "Long axis orientation (degrees)", 
    ylabel = "fraction", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)

smallDP = CairoMakie.stairs!(ax18, sort(filtered_df_small[!, :MeanOrientationCtoDP]), ecdf(sort(filtered_df_small[!, :MeanOrientationCtoDP]))(sort(filtered_df_small[!, :MeanOrientationCtoDP])), color = (Makie.wong_colors()[1], 0.5))
bigDP = CairoMakie.stairs!(ax18, sort(filtered_df_big[!, :MeanOrientationCtoDP]), ecdf(sort(filtered_df_big[!, :MeanOrientationCtoDP]))(sort(filtered_df_big[!, :MeanOrientationCtoDP])), color = (Makie.wong_colors()[1], 1))
smallRachis = CairoMakie.stairs!(ax18, sort(filtered_df_small[!, :MeanAnglesLongAxisRachis]), ecdf(sort(filtered_df_small[!, :MeanAnglesLongAxisRachis]))(sort(filtered_df_small[!, :MeanAnglesLongAxisRachis])), color = (Makie.wong_colors()[3], 0.5))
bigRachis = CairoMakie.stairs!(ax18, sort(filtered_df_big[!, :MeanAnglesLongAxisRachis]), ecdf(sort(filtered_df_big[!, :MeanAnglesLongAxisRachis]))(sort(filtered_df_big[!, :MeanAnglesLongAxisRachis])), color = (Makie.wong_colors()[3], 1))
random = CairoMakie.lines!(ax18, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1,  color = Makie.wong_colors()[2])

PolarHist(filtered_df_small[!, :MeanOrientationCtoDP], f18, "Long axis\n to DP axis, small cells", 1, 0.3, 3)
PolarHist(filtered_df_big[!, :MeanOrientationCtoDP], f18, "Long axis \n to DP axis, big cells", 1, 0.75, 4)
PolarHist(filtered_df_small[!, :MeanAnglesLongAxisRachis], f18, "Long axis\n to rachis normal, small cells", 3, 0.3, 5)
PolarHist(filtered_df_big[!, :MeanAnglesLongAxisRachis], f18, "Long axis \n to rachis normal, big cells", 3, 0.75, 6)

CairoMakie.Legend(f18[1,2], [smallDP, bigDP, smallRachis, bigRachis, random], ["cell volume < 125 µm³, DP", "cell volume ≥ 125 µm³, DP", "cell volume < 125 µm³, Rachis", "cell volume ≥ 125 µm³, Rachis","theoretical random"])
CairoMakie.display(f18)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/LongAxisOrientionToDPandRachisSmallVsBigCells.pdf", f18, px_per_unit = 2)

f20 = CairoMakie.Figure(size = (1350, 450))
ax20 = CairoMakie.Axis(f20[1,1], 
    title = " Cells long axis orientation \n at mitosis, n = $(nrow(Mitosis_df))", 
    xlabel = "Cells long axis orientation (degrees)", 
    ylabel = "fraction", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)

DP = CairoMakie.stairs!(ax20, sort(Mitosis_df[!, :MeanOrientationCtoDP]), ecdf(sort(Mitosis_df[!, :MeanOrientationCtoDP]))(sort(Mitosis_df[!, :MeanOrientationCtoDP])), color = (Makie.wong_colors()[1], 1))
rachis = CairoMakie.stairs!(ax20, sort(Mitosis_df[!, :MeanAngleLongAxisRachis]), ecdf(sort(Mitosis_df[!, :MeanAngleLongAxisRachis]))(sort(Mitosis_df[!, :MeanAngleLongAxisRachis])), color = (Makie.wong_colors()[3], 1))
random = CairoMakie.lines!(ax20, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1, color = (Makie.wong_colors()[2], 1))    
CairoMakie.Legend(f20[1,2], [DP, rachis, random], ["D/P axis", "rachis surface normal", "theoretical random"])

PolarHist(Mitosis_df[!, :MeanOrientationCtoDP], f20, "Long axis axis \n to DP axis", 1, 0.75, 3)
PolarHist(Mitosis_df[!, :MeanAngleLongAxisRachis], f20, "Long axis axis \n to rachis normal", 3, 0.75, 4)

CairoMakie.display(f20)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/CellsLongAxisOrientationMitosis.pdf", f20) 


f21 = CairoMakie.Figure(size = (2250, 450))
ax21 = CairoMakie.Axis(f21[1,1], 
    title = " Cells long axis orientation to D/P axis and rachis normal \n at interphase, n (round) = $(nrow(Interphase_small_prolate)), n (elongated) = $(nrow(Interphase_big_prolate))",  
    xlabel = "Cells long axis orientation (degrees)", 
    ylabel = "fraction", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)
DPaxisSmall = CairoMakie.stairs!(ax21, sort(Interphase_small_prolate[!, :MeanOrientationCtoDP]), ecdf(sort(Interphase_small_prolate[!, :MeanOrientationCtoDP]))(sort(Interphase_small_prolate[!, :MeanOrientationCtoDP])), color = (Makie.wong_colors()[1], 0.5))
DPaxisBig = CairoMakie.stairs!(ax21, sort(Interphase_big_prolate[!, :MeanOrientationCtoDP]), ecdf(sort(Interphase_big_prolate[!, :MeanOrientationCtoDP]))(sort(Interphase_big_prolate[!, :MeanOrientationCtoDP])), color = (Makie.wong_colors()[1], 1))

rachisSmall = CairoMakie.stairs!(ax21, sort(Interphase_small_prolate[!, :MeanAnglesLongAxisRachis]), ecdf(sort(Interphase_small_prolate[!, :MeanAnglesLongAxisRachis]))(sort(Interphase_small_prolate[!, :MeanAnglesLongAxisRachis])), color = (Makie.wong_colors()[3], 0.5))
rachisBig = CairoMakie.stairs!(ax21, sort(Interphase_big_prolate[!, :MeanAnglesLongAxisRachis]), ecdf(sort(Interphase_big_prolate[!, :MeanAnglesLongAxisRachis]))(sort(Interphase_big_prolate[!, :MeanAnglesLongAxisRachis])), color = (Makie.wong_colors()[3], 1))

random = CairoMakie.lines!(ax21, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1, color = (Makie.wong_colors()[2], 1))

PolarHist(Interphase_small_prolate[!, :MeanOrientationCtoDP], f21, "Long axis axis \n to DP axis, round cells", 1, 0.3, 3)
PolarHist(Interphase_big_prolate[!, :MeanOrientationCtoDP], f21, "Long axis axis \n to DP axis, elongated cells", 1, 0.75, 4)
PolarHist(Interphase_small_prolate[!, :MeanAnglesLongAxisRachis], f21, "Long axis axis \n to rachis normal, round cells", 3, 0.3, 5)
PolarHist(Interphase_big_prolate[!, :MeanAnglesLongAxisRachis], f21, "Long axis axis \n to rachis normal, elongated", 3, 0.75, 6)

CairoMakie.Legend(f21[1,2], [DPaxisSmall, DPaxisBig,rachisSmall,  rachisBig, random], ["D/P axis\n prolate index < 0.4", "D/P axis\n prolate index ≥ 0.4", "rachis surface normal \n prolate index < 0.4", "rachis surface normal \n prolate index ≥ 0.4", "theoretical random"])
CairoMakie.display(f21)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/CellsLongAxisOrientationInterphaseBigvsSmallProlate.pdf", f21, px_per_unit = 2) 

f22 = CairoMakie.Figure(size = (700, 450))
ax22 = CairoMakie.Axis(f22[1,1], 
    title = "Interphase cells prolate index according to cell volume,\n algorithm = PseudorandomJitter, \n n (small) = $(nrow(filtered_df_small)), n (big) = $(nrow(filtered_df_big))",  
    xlabel = "Cell volume", 
    ylabel = "Prolate index", 
    xticks = ([1,2], ["< 125 µm³", "≥ 125 µm³"]), 
    xgridvisible  = false,
    ygridvisible  = false,
    limits = (0.5, 2.5, 0, 1),
    )
Small = SwarmMakie.beeswarm!(ax22, (filtered_df_small[!, :CellSize] .== "Big").+1, filtered_df_small[!, :MeanProlateIndex], algorithm = UniformJitter( jitter_width = 40), color = :blue, alpha = 0.5, markersize = 4)
Big = SwarmMakie.beeswarm!(ax22,(filtered_df_big[!, :CellSize] .== "Big").+1, filtered_df_big[!, :MeanProlateIndex], algorithm = UniformJitter( jitter_width = 40), color = :green, alpha = 0.5, markersize = 4)
CairoMakie.Legend(f22[1,2], [Small, Big], ["cell volume < 125 µm³", "cell volume ≥ 125 µm³"])
CairoMakie.display(f22)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/InterphaseCellsProlateIndexPseudorandomJitter.pdf", f22, px_per_unit = 2) 


MeanInterphasePA =  mean(Interphase_df[sample_rows, :MeanProlateIndex])
StdInterphasePA =  std(Interphase_df[sample_rows, :MeanProlateIndex])
MeanMitosisPA =  mean(Mitosis_df[!, :MeanProlateIndex])
StdMitosisPA =  std(Mitosis_df[!, :MeanProlateIndex])



f23 = CairoMakie.Figure(size = (500, 600))
ax23 = CairoMakie.Axis(f23[1,1], 
    title = "cells prolate index, \n n (Interphase) = $(nrow(Interphase_df)), n (mitosis) = $(nrow(Mitosis_df))",  
    ylabel = "Prolate index", 
    xticks = ([1,2], ["Interphase", "Mitosis"]), 
    yticks = collect(0:0.2:1),
    xgridvisible  = false,
    ygridvisible  = false,
    limits = (0.5, 2.5, 0, 1),
    )
SwarmMakie.beeswarm!(ax23, fill(1, nrow(Interphase_df)), Interphase_df[!, :MeanProlateIndex], 
color = (Makie.wong_colors()[5],0.6),
markersize = 10,
)
SwarmMakie.beeswarm!(ax23, fill(2, nrow(Mitosis_df)), Mitosis_df[!, :MeanProlateIndex], 
color = (Makie.wong_colors()[6], 0.5),
#markersize = 10
)

errorbars!( ax23, [1.0], [MeanInterphasePA], [StdInterphasePA], [StdInterphasePA], color = :black, whiskerwidth = 10)
errorbars!(ax23, [2.0], [MeanMitosisPA], [StdMitosisPA], [StdMitosisPA], color = :black, whiskerwidth = 10)
crossbar!( ax23, [1.0], [MeanInterphasePA], [0], [0])
crossbar!( ax23, [2.0], [MeanMitosisPA], [0], [0])

CairoMakie.display(f23)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/ProlateIndexRainCloudNotSeparared.pdf", f23, px_per_unit = 20) 

CellShapeCleanedSizeMitosis = KeepBigCells(Interphase_df, Mitosis_df)
Mitosis_df = dropmissing(Mitosis_df, [:MeanCellVolume])
sample_rows = sample(1:nrow(CellShapeCleanedSizeMitosis), 102, replace=false)
TTest = HypothesisTests.UnequalVarianceTTest(CellShapeCleanedSizeMitosis[sample_rows, :MeanProlateIndex], Mitosis_df[!, :MeanProlateIndex])
TTestV = HypothesisTests.UnequalVarianceTTest(CellShapeCleanedSizeMitosis[sample_rows, :MeanCellVolume], Mitosis_df[!, :MeanCellVolume])

Interphase102Random = copy(CellShapeCleanedSizeMitosis[sample_rows, [:Gonad, :TrackID, :Duration, :MeanCellVolume, :MeanProlateIndex, :MeanDistanceDTC,
                                                                :MeanOrientationCtoDP, :MeanAnglesLongAxisRachis,
                                                                :MeanCellEllipsoidAxisCX, :MeanCellEllipsoidAxisCY, :MeanCellEllipsoidAxisCZ,
                                                                :MeanCellEllipsoidAxisBX, :MeanCellEllipsoidAxisBY, :MeanCellEllipsoidAxisBZ,
                                                                :MeanCellEllipsoidAxisAX, :MeanCellEllipsoidAxisAY, :MeanCellEllipsoidAxisAZ,
                                                                :MeanAngleProjectLongestAxisDP,
                                                                ]
                                                 ])


rename!(Interphase102Random, [:Gonad, :TrackID, :Duration, :CellVolume, :ProlateIndex, :DistanceDTC,
                              :OrientationLongAxisToDP, :OrientationLongAxisToRachis, 
                              :LongAxisX, :LongAxisY, :LongAxisZ,
                              :MediumAxisX, :MediumAxisY, :MediumAxisZ,
                              :ShortAxisX, :ShortAxisY, :ShortAxisZ,
                              :OrientationLongestProjectedAxisToDP
                              ])

CSV.write("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/Interphase102Random.csv", Interphase102Random, delim = ",")

MeanInterphase =  mean(CellShapeCleanedSizeMitosis[sample_rows, :MeanProlateIndex])
StdInterphase =  std(CellShapeCleanedSizeMitosis[sample_rows, :MeanProlateIndex])
MeanMitosis =  mean(Mitosis_df[!, :MeanProlateIndex])
StdMitosis =  std(Mitosis_df[!, :MeanProlateIndex])

f24 = CairoMakie.Figure(size = (600, 600))
ax24 = CairoMakie.Axis(f24[1,1], 
    title = "Interphase and mitotic cells prolate index, \n n (Interphase all) = $(nrow(Interphase_df)), n (Interphase) = $(nrow(CellShapeCleanedSizeMitosis[sample_rows,:])), n (Mitotic) = $(nrow(Mitosis_df)) \n p value =  $(pvalue(TTest))\n mean Interphase = $(MeanInterphase) +/- $(StdInterphase)\n mean Mitosis = $(MeanMitosis) +/- $(StdMitosis)",  
    xlabel = "Cells", 
    ylabel = "Prolate index", 
    xticks = ([1,2,3], ["Interphase (all)", "Interphase", "Mitosis"]),
    yticks = collect(0:0.2:1),
    xgridvisible  = false,
    ygridvisible  = false,
    limits = (0.5, 3.5, 0, 1),
    )

InterphaseAll = violin!(ax24,  fill(1, nrow(Interphase_df)), Interphase_df[!, :MeanProlateIndex], color = (Makie.wong_colors()[5],0.5))
violin!(ax24,  fill(2, nrow(CellShapeCleanedSizeMitosis[sample_rows,:])), CellShapeCleanedSizeMitosis[sample_rows, :MeanProlateIndex], color =  (Makie.wong_colors()[5], 0.25))
Interphase = SwarmMakie.beeswarm!(ax24,  fill(2, nrow(CellShapeCleanedSizeMitosis[sample_rows,:])), CellShapeCleanedSizeMitosis[sample_rows, :MeanProlateIndex], algorithm = SimpleBeeswarm(), color =  Makie.wong_colors()[5], alpha = 0.75)
violin!(ax24,fill(3, nrow(Mitosis_df)), Mitosis_df[!, :MeanProlateIndex], color =  (Makie.wong_colors()[6], 0.25))

Mitosis = SwarmMakie.beeswarm!(ax24,fill(3, nrow(Mitosis_df)), Mitosis_df[!, :MeanProlateIndex], algorithm = SimpleBeeswarm(), color =  Makie.wong_colors()[6], alpha = 0.75)
errorbars!( ax24, [2.0], [MeanInterphase], [StdInterphase], [StdInterphase], color = :black, whiskerwidth = 10)
errorbars!(ax24, [3.0], [MeanMitosis], [StdMitosis], [StdMitosis], color = :black, whiskerwidth = 10)
crossbar!( ax24, [2.0], [MeanInterphase], [0], [0])
crossbar!( ax24, [3.0], [MeanMitosis], [0], [0])
CairoMakie.Legend(f24[1,2], [Interphase, Mitosis], ["Interphase cells", "Mitotic cells"])
CairoMakie.display(f24)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/ProlateIndexMitosisVsInterphaseAllAndSelected.pdf", f24, px_per_unit = 2) 


MeanInterphaseV =  mean(CellShapeCleanedSizeMitosis[sample_rows, :MeanCellVolume])
StdInterphaseV =  std(CellShapeCleanedSizeMitosis[sample_rows, :MeanCellVolume])
MeanMitosisV =  mean(Mitosis_df[!, :MeanCellVolume])
StdMitosisV =  std(Mitosis_df[!, :MeanCellVolume])


f25 = CairoMakie.Figure(size = (600, 600))
#f25 = CairoMakie.Figure()
ax25 = CairoMakie.Axis(f25[1,1], 
    title = "Interphase and mitotic cells volume, \n n (Interphase all) = $(nrow(Interphase_df)), n (Interphase) = $(nrow(CellShapeCleanedSizeMitosis[sample_rows,:])), n (Mitotic) = $(nrow(Mitosis_df)) \n p value =  $(pvalue(TTestV))\n mean Interphase = $(MeanInterphaseV) +/- $(StdInterphaseV)\n mean Mitosis = $(MeanMitosisV) +/- $(StdMitosisV)",  
    xlabel = "Cells", 
    ylabel = "Cell volume (µm³)", 
    xticks = ([1,2,3], ["Interphase (all)", "Interphase", "Mitosis"]), 
    yticks = collect(50:25:250),
    xgridvisible  = false,
    ygridvisible  = false,
    limits = (0.5, 3.5, 50, 250),
    )
InterphaseAll = violin!(ax25,  fill(1, nrow(Interphase_df)), Interphase_df[!, :MeanCellVolume], color = (Makie.wong_colors()[5],0.5))
violin!(ax25,  fill(2, nrow(CellShapeCleanedSizeMitosis[sample_rows,:])), CellShapeCleanedSizeMitosis[sample_rows, :MeanCellVolume], color =  (Makie.wong_colors()[5], 0.25))
Interphase = SwarmMakie.beeswarm!(ax25,  fill(2, nrow(CellShapeCleanedSizeMitosis[sample_rows,:])), CellShapeCleanedSizeMitosis[sample_rows, :MeanCellVolume], algorithm = SimpleBeeswarm(), color =  Makie.wong_colors()[5], alpha = 0.75)
violin!(ax25,fill(3, nrow(Mitosis_df)), Mitosis_df[!, :MeanCellVolume], color =  (Makie.wong_colors()[6], 0.25))

Mitosis = SwarmMakie.beeswarm!(ax25,fill(3, nrow(Mitosis_df)), Mitosis_df[!, :MeanCellVolume], algorithm = SimpleBeeswarm(), color =  Makie.wong_colors()[6], alpha = 0.75)
errorbars!( ax25, [2.0], [MeanInterphaseV], [StdInterphaseV], [StdInterphaseV], color = :black, whiskerwidth = 10)
errorbars!(ax25, [3.0], [MeanMitosisV], [StdMitosisV], [StdMitosisV], color = :black, whiskerwidth = 10)
crossbar!( ax25, [2.0], [MeanInterphaseV], [0], [0])
crossbar!( ax25, [3.0], [MeanMitosisV], [0], [0])
CairoMakie.Legend(f25[1,2], [Interphase, Mitosis], ["Interphase cells", "Mitotic cells"])
CairoMakie.display(f25)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/CellVolumeMitosisVsInterphaseAllAndSelected.pdf", f25, px_per_unit = 2) 


f26 = CairoMakie.Figure(size = (600, 600))
ax26 = CairoMakie.Axis(f26[1,1], 
    title = "Interphase and mitotic cells volume, \n n (Interphase) = $(nrow(Interphase_df)), n (Mitotic) = $(nrow(Mitosis_df))",  
    xlabel = "Cells", 
    ylabel = "Cell volume (µm³)", 
    xticks = ([1,2], ["Interphase", "Mitosis"]), 
    yticks = collect(0:50:300),
    xgridvisible  = false,
    ygridvisible  = false,
    limits = (0.5, 2.5, 0, 300),
    )
Interphase = rainclouds!(ax26,  fill(1, nrow(Interphase_df)), Interphase_df[!, :MeanCellVolume], color =  Makie.wong_colors()[5])
Mitosis = rainclouds!(ax26,fill(2, nrow(Mitosis_df)), Mitosis_df[!, :MeanCellVolume], color =  Makie.wong_colors()[6])
CairoMakie.display(f26)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/CellVolumeMitosisVsInterphaseSelectAll.pdf", f26, px_per_unit = 2) 

f27 = CairoMakie.Figure(size = (1350, 450));
ax27 = CairoMakie.Axis(f27[1,1], 
    title = "Longest projected cell axis on rachis surface \n orientation to projected DP axis \n mitosis n = $(nrow(Mitosis_df)), interphase n =$(nrow(Interphase_df))", 
    xlabel = "Orientation to projected DP axis(degrees)", 
    ylabel = "Fraction", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)

Mitosis = CairoMakie.stairs!(ax27, sort(Mitosis_df[!, :MeanAngleProjectLongestAxisDP]), ecdf(sort(Mitosis_df[!, :MeanAngleProjectLongestAxisDP]))(sort(Mitosis_df[!, :MeanAngleProjectLongestAxisDP])), color = Makie.wong_colors()[6])
Interphase = CairoMakie.stairs!(ax27, sort(Interphase_df[!, :MeanAngleProjectLongestAxisDP]), ecdf(sort(Interphase_df[!, :MeanAngleProjectLongestAxisDP]))(sort(Interphase_df[!, :MeanAngleProjectLongestAxisDP])), color = Makie.wong_colors()[5])
random = CairoMakie.lines!(ax27, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1, color = Makie.wong_colors()[7])    
CairoMakie.Legend(f27[1,2], [Mitosis, Interphase, random], ["Mitosis","Interphase", "theoretical random"])
PolarHist(Mitosis_df[!, :MeanAngleProjectLongestAxisDP], f27, "Mitosis", 6, 0.75, 3)
PolarHist(Interphase_df[!, :MeanAngleProjectLongestAxisDP], f27, "Interphase", 5, 0.75, 4) 
CairoMakie.display(f27)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/CellsProjectedLongestAxisOrientationDPaxiswMitosisandInterphase.pdf", f27, px_per_unit = 2) 

f28 = CairoMakie.Figure(size = (1100, 450));
ax28 = CairoMakie.Axis(f28[1,1], 
    title = "Projected mitotic spindle on rachis surface \n orientation to longest projected axis \n  n = $(nrow(Mitosis_df))",
    ylabel = "Fraction", 
    xlabel = "Orientation to longest projected axis (degrees)", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)

Mitosis = CairoMakie.stairs!(ax28, sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis]), ecdf(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis]))(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis])), color = Makie.wong_colors()[6])
random = CairoMakie.lines!(ax28, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1, color = Makie.wong_colors()[7])    
CairoMakie.Legend(f28[1,2], [Mitosis, random], ["Projected spindle at AO", "theoretical random"])
PolarHist(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis], f28, "Projected Spindle", 6, 0.75, 3)
CairoMakie.display(f28)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/ProjectedSpindleToLongestProjectedAxis.pdf", f28, px_per_unit = 2) 

f29 = CairoMakie.Figure(size = (1100, 450));
ax29 = CairoMakie.Axis(f29[1,1], 
    title = "Projected mitotic spindle on rachis surface \n orientation to shortest projected axis \n  n = $(nrow(Mitosis_df))",
    ylabel = "Fraction", 
    xlabel = "Orientation to shortest projected axis (degrees)", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)

Mitosis = CairoMakie.stairs!(ax29, sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectShortAxis]), ecdf(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis]))(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis])), color = Makie.wong_colors()[6])
random = CairoMakie.lines!(ax29, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1, color = Makie.wong_colors()[7])    
CairoMakie.Legend(f29[1,2], [Mitosis, random], ["Projected spindle at AO", "theoretical random"])
PolarHist(Mitosis_df[!, :MeanAngleProjectSpindleProjectShortAxis], f29, "Projected Spindle", 6, 0.75, 3)
CairoMakie.display(f29)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/ProjectedSpindleToShortestProjectedAxis.pdf", f29, px_per_unit = 2) 


f30 = CairoMakie.Figure(size = (650, 500));
ax30 = CairoMakie.Axis(f30[1,1], 
    title = "Mitotic spindle orientation \n  n = $(nrow(Mitosis_df))",
    xlabel = "Rachis surface normal (degrees)", 
    ylabel = "Cell long axis (degrees)", 
    #zlabel = "Orientation to DP axis (degrees)",
    limits = (60, 90, 0, 90),
    xgridvisible  = false,
    ygridvisible  = false)

Mitosis =  GLMakie.scatter!(ax30, 
Mitosis_df[!, :MeanAngleSpindleRachisAnaphaseOnset], 
Mitosis_df[!, :MeanAngleSpindleLongAxisAnaphaseOnset],
#Mitosis_df[!, :MeanAngleSpindleDPaxis]
)
CairoMakie.display(f30)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/MitoticSpindleOrientationRachisvsLongAxis.pdf", f30, px_per_unit = 2) 

f31 = CairoMakie.Figure(size = (650, 500));
ax31 = CairoMakie.Axis(f31[1,1], 
    title = "Mitotic spindle orientation \n  n = $(nrow(Mitosis_df))",
    xlabel = "Orientation to long axis (degrees)", 
    ylabel = "Orientation to D/P axis (degrees)", 
    #zlabel = "Orientation to DP axis (degrees)",
    limits = (0, 90, 0, 90),
    xgridvisible  = false,
    ygridvisible  = false)

Mitosis =  GLMakie.scatter!(ax31, 
Mitosis_df[!, :MeanAngleSpindleLongAxisAnaphaseOnset], 
Mitosis_df[!, :MeanAngleSpindleDPaxis],
#Mitosis_df[!, :MeanAngleSpindleDPaxis]
)
CairoMakie.display(f31)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/MitoticSpindleOrientationLongAxisvsDPaxis.pdf", f31, px_per_unit = 2) 


f32 = CairoMakie.Figure(size = (1350, 450));
ax32 = CairoMakie.Axis(f32[1,1], 
    title = "Projected mitotic spindle on rachis surface \n orientation to projected DP axis and longest projected axis \n  n = $(nrow(Mitosis_df))",
    ylabel = "Fraction", 
    xlabel = "Projected mitotic Spindle orientation (degrees)", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)

DPAxis = CairoMakie.stairs!(ax32, sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectDPaxis]), ecdf(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectDPaxis]))(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectDPaxis])), color = Makie.wong_colors()[6])
LongAxis = CairoMakie.stairs!(ax32, sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis]), ecdf(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis]))(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis])), color = Makie.wong_colors()[1])

random = CairoMakie.lines!(ax32, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1, color = Makie.wong_colors()[7])    
CairoMakie.Legend(f32[1,2], [DPAxis, LongAxis, random], ["Projected DP axis", "Longest projected axis", "theoretical random"])
PolarHist(Mitosis_df[!, :MeanAngleProjectSpindleProjectDPaxis], f32, "Projected D/P axis", 6, 0.75, 3)
PolarHist(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis], f32, "Longest projected axis", 1, 0.75, 4)

CairoMakie.display(f32)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/ProjectedSpindleToProjetedDPaxisAndProjectedLongestAxis.pdf", f32, px_per_unit = 2) 



f33 = CairoMakie.Figure(size = (1100, 450));
ax33 = CairoMakie.Axis(f33[1,1], 
    title = "Longest projected axis \n orientation to projected DP axis \n  n = $(nrow(Mitosis_df))",
    ylabel = "Fraction", 
    xlabel = "Orientation to longest projected axis (degrees)", 
    limits = (0, 90, 0, 1),
    xticks = 0:10:90,
    xgridvisible  = false,
    ygridvisible  = false,
    yticks = 0:0.1:1)

Mitosis = CairoMakie.stairs!(ax33, sort(Mitosis_df[!, :MeanAngleProjectDPaxisProjectLongAxis]), ecdf(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis]))(sort(Mitosis_df[!, :MeanAngleProjectSpindleProjectLongAxis])), color = Makie.wong_colors()[6])
random = CairoMakie.lines!(ax33, collect(0:0.5:90), -cos.(collect(0:0.5:90).*pi/180) .+1, color = Makie.wong_colors()[7])    
CairoMakie.Legend(f33[1,2], [Mitosis, random], ["Longest projected axis to projected DP axis", "theoretical random"])
PolarHist(Mitosis_df[!, :MeanAngleProjectDPaxisProjectLongAxis], f33, "Longest projected axis \n to projected DP axis", 6, 0.75, 3)
CairoMakie.display(f33)
save("M:/Labbe/VincentPoupart/Papers/2024/RedaSpindleOrientation/JuliaFigures/ProjectedLongesttoProjectedAxisDPaxis.pdf", f33, px_per_unit = 2) 