

function align_vectors(vectors, grid_length) 
    # takes as inpput a list of variable length vectors 
    # and aligns then to a grid with 1:3 intrepolation

    # Determine the range of the new uniform grid
    grid = collect(1:grid_length)
    
    # Align all vectors to the new grid
    aligned_vectors = [interpolate_vector(vec, grid) for vec in vectors]

    return hcat(aligned_vectors...)
end

# Linear interpolation function
function interp1d(x, y,xi)
    if xi <= x[1]
        return y[1]
    elseif xi >= x[end]
        return y[end]
    else
        for i in 1:length(x)-1
            if x[i] <= xi <= x[i+1]
                return y[i] + (y[i+1] - y[i]) * (xi - x[i]) / (x[i+1] - x[i])
            end
        end
    end
end

function interpolate_vector(vector, grid)
    # Create the original grid for the vector
    original_grid =  collect(1:length(vector))

    # Perform linear interpolation
    interp_vector = [interp1d(original_grid, vector, x) for x in grid]
    return interp_vector
end
