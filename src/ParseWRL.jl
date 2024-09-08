function skip_pass(marker, lines)
    # Skip until reaching the line which contains the marker, then stop
    for line in lines
        if occursin(marker, line)
            # Include the line containing the marker and set include_line to true
            #println(line)
            return IterTools.chain([line], lines)
        end
    end
    # If marker not found, return all lines
    return lines
end

function take(lines)
    #=
    Take and return those lines which contain a marker
    =#
    result = []
    for line in lines
        if occursin("]", line) == false
            push!(result, line)
        else
            push!(result, line)
            println(line)
            break
        end
    end
    return result
end

#=
function take(marker, lines)
    #=
    Take and return those lines which contain a marker
    =#
    result = []
    for line in lines
        if occursin(marker, line) 
            push!(result, line)
        else
            break
        end
    end
    return result
end
=#
function CatData(fazers)
    type = Vector{Vector{Float64}}()
    for line in fazers
        #matches = eachmatch(r"[-+]?\d*\.?\d+", line)
        matches = eachmatch(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", line)
        parsed = [parse(Float64, m.match) for m in matches]
        push!(type, parsed)
    end
    return type
end

function CatDataIndices(fazers)
    type = Vector{Tuple}()
    for line in fazers
        tokens = filter(x -> x != "-1", split(strip(line, ['[', ']']), r"\s*[,[:alpha:]\[\]\s]+"))
        tokens = filter(x -> x != "", tokens)
        indices = Vector{Int64}()
        for token in tokens
            idx = parse(Int64, token) + 1 # add 1 to parsed value because Julia indices start at 1, not 0
            push!(indices, idx)
        end
        # Splitting coord_lines into groups of 3 indices
        for idx in 1:3:length(indices)-2
            push!(type, Tuple(indices[idx:idx+2]))
        end
    end
    return type
end




function CatDataFacets(coord_index, verts)
    facets = Vector{Vector{Vector{Float64}}}()
    centroids = Vector{Vector{Float64}}()
    areas = Vector{Float64}()
    normals = Vector{Vector{Float64}}()
    for indices in coord_index
        face = Vector{Vector{Float64}}()
        for i in indices
            push!(face, verts[i])
        end
        
        push!(facets, face)
     
        try
        push!(centroids, mean(face))
        catch err
            println("indices = ", indices)
            println(verts[indices[2]])
        end
        VA = face[2] - face[1]
        VB = face[3] - face[1]
        push!(areas, 0.5(norm(cross(VA, VB))))
        push!(normals, (cross(VA, VB)))
    end
    return facets, centroids, areas, normals
end

function parse_indexed_face_set(lines)
    ###### Parse one block of 'geometry IndexedFaceSet'
    # Parse the vertices
    lines1 = skip_pass("point", lines)
    #point_lines = take(" ", lines1)
    point_lines = take(lines1)
    verts = CatData(point_lines)
    # Parse the normal vectors to each vertex 
    # Reset lines iterator
    lines3 = skip_pass("vector", lines)
    vector_lines = take(lines3)
    vectors = CatData(vector_lines)
    # Parse the coordinate indices of the facets
    # Reset lines iterator
    lines2 = skip_pass("coordIndex", lines)
    coord_lines = take(lines2)
    coord_index = CatDataIndices(coord_lines)
    # Using coord_indices to build facets
    facets, centroids, areas, normals = CatDataFacets(coord_index, verts)
    set = Dict("verts" => verts, "facets" => facets, "coordIndex" => coord_index,"vectors" => vectors, "normals" => normals, "centroids" => centroids, "areas" => areas)
    return set
end


function ExtractDataVRML(WRLpath)
    IndexedFaceSets = Vector{Dict{String, Vector{Any}}}()
    open(WRLpath) do infile
        i = 0
        for line in eachline(infile)
            if occursin("IndexedFaceSet", line)
                set = parse_indexed_face_set(eachline(infile))
                push!(IndexedFaceSets, set)
                set = nothing
                i += 1
                println("timepoint $i is done")
            end
        end
    end
    return IndexedFaceSets
end
