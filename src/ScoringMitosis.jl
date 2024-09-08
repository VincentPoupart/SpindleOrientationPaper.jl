
function PlotSplindleLengthVsFrame(Celloutput::Dict{String,DataFrame}, cellule::String)
    df = Celloutput[cellule]
    frame = df[!, :FRAME]
    spindle_length = df[!, :Spindle_Length]
    f = Figure(size=(800, 600))
    ax1 = Axis(f[2, 1], title="$(cellule)", limits=(nothing, (0, 12)))
    lines!(ax1, frame, spindle_length, color=:blue, linewidth=3)
    scatter!(ax1, frame, spindle_length, markersize=30, marker=:circle, color=:transparent, strokewidth=3, strokecolor=:blue)
    ax = content(f[2, 1])
    Makie.deactivate_interaction!(ax, :rectanglezoom)
    spoint = select_point(ax.scene, marker=:circle)
    done = Button(f[3, 1]; label="done", tellwidth=false)
    redo = Button(f[1, 1]; label="redo", tellwidth=false)
    isdone = Observable(false)
    isredo = Observable(false)
    on(done.clicks) do clicks
        isdone[] = true
    end
    on(redo.clicks) do clicks
        isredo[] = true
    end
    MousePress = []
    NumberPress = Observable(0)
    NEBD = Observable(NaN)
    StartCongression = Observable(NaN)
    AnaphaseOnset = Observable(NaN)
    on(events(ax).mousebutton) do event
        if event.action == Mouse.press
            push!(MousePress, true)
            NumberPress[] = length(MousePress)
            if event.button == Mouse.left
                on(spoint) do z
                    x, y = z
                    xb = findnearest(frame, x)
                    if frame[xb] < x
                        xc = frame[xb+1]
                    else
                        xc = frame[xb]
                    end
                    if NumberPress[] == 1
                        MitoticStage = "NEBD"
                        couleur = :blue
                        NEBD[] = xc
                    elseif NumberPress[] == 2
                        MitoticStage = "Start of Congression"
                        couleur = :red
                        StartCongression[] = xc
                    elseif NumberPress[] == 3
                        MitoticStage = "Anaphase Onset"
                        couleur = :green
                        AnaphaseOnset[] = xc
                    end
                    if NumberPress[] <= 3
                        vlines!(ax1, xc; linewidth=3, linestyle=:dashdot, color=couleur)
                        text!(ax1, xc - 0.25, 6, text=MitoticStage; rotation=pi / 2, color=couleur, fontsize=30)
                    end
                end
            end
        end
    end
    return f, spoint, isdone, isredo, NumberPress, NEBD, StartCongression, AnaphaseOnset
end


#Running throughout the celloutput and scoring mitosis for each cells
function ClickStepsMitosis(Celloutput::Dict)
    GLMakie.activate!(; visible=true, focus_on_show=true, fullscreen = true, render_on_demand = true)
    Scoring = Dict{String,DataFrame}()
    k = collect(keys(Celloutput))
    for i in eachindex(k)
        f, spoint, isdone, isredo, NumberPress, NEBD, StartCongression, AnaphaseOnset = PlotSplindleLengthVsFrame(Celloutput, k[i])
        display(f, scalefactor = 2)
        while isdone[] == false
            sleep(0.5)
            isopen(f.scene) || break
            if isredo[] == true || NumberPress[] > 3
                f, spoint, isdone, isredo, NumberPress, NEBD, StartCongression, AnaphaseOnset = PlotSplindleLengthVsFrame(Celloutput, k[i])
                display(f, scalefactor = 2)
            end
        end
        isopen(f.scene) || break
        df = DataFrame(NEBD=NEBD[], StartCongression=StartCongression[], AnaphaseOnset=AnaphaseOnset[])
        Scoring[k[i]] = df
    end
    Scoring = sort(Scoring)
    k = collect(keys(Scoring))
    scoring_df = DataFrame()
    Cell = Vector{Union{String,Missing}}(undef, length(Scoring))
    NEBD = Vector{Float64}(undef, length(Scoring))
    StartCongression = Vector{Float64}(undef, length(Scoring))
    AnaphaseOnset = Vector{Float64}(undef, length(Scoring))
    for i in eachindex(k)
        Cell[i] = k[i]
        NEBD[i] = Scoring[k[i]][1, :NEBD]
        StartCongression[i] = Scoring[k[i]][1, :StartCongression]
        AnaphaseOnset[i] = Scoring[k[i]][1, :AnaphaseOnset]
    end
    scoring_df[!, :Cell] = Cell
    scoring_df[!, :NEBD] = NEBD
    scoring_df[!, :StartCongression] = StartCongression
    scoring_df[!, :AnaphaseOnset] = AnaphaseOnset
    return scoring_df
end
