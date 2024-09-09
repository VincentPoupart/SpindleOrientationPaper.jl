module SpindleOrientationPaper


using DataFrames
using CSV
using DataFramesMeta
using StatsBase
using LinearAlgebra
using Statistics
using CairoMakie
using GLMakie
using HypothesisTests
using SwarmMakie
using PlotlyJS



export
    SpindleOrientationCellShape,
    ClickStepsMitosis,
    ExtractDataVRML,
    GeneratingFigures

include("AddNormalVectors.jl")
include("AddProjectAxisCellShape.jl")
include("CellShape.jl")
include("CellShapeTracksFinder.jl")
include("DataFetching.jl")
include("GeneratingFigures.jl")
include("ImportSavedData.jl")
include("ParseWRL.jl")
include("PlotSurface.jl")
include("ScoringMitosis.jl")
include("SmallFunctions.jl")

end
