"""
    create_velocity_dis(parameters::parameters) -> Distribution

Creates a velocity distribution for agents based on the simulation parameters.

## Arguments
- `parameters::parameters`: A struct containing the simulation parameters, which must include:
  - `velocity_dis`: The type of distribution to use (`Uniform` or `Normal`).
  - `velocity_variance`: The variance of the velocity distribution.
  - `velocity`: The mean velocity of the agents.

## Returns
- `Distribution`: A truncated distribution object representing the velocity distribution of agents.
     The distribution is truncated at the lower bound of 0 to ensure that velocity values are non-negative.

## Behavior
- If both `velocity` and `velocity_variance` are 0, a `Normal` distribution with mean 0 and variance 0 is returned.
- If `velocity_dis` is `Uniform`, a truncated `Uniform` distribution ranging from `(velocity - variance)` to `(velocity + variance)` is created.
- If `velocity_dis` is `Normal`, a truncated `Normal` distribution with the specified mean (`velocity`) and variance (`velocity_variance`) is created.
- The truncation ensures that the resulting velocity values are always non-negative, as negative velocities are not physically meaningful in this context.

## Example
```julia
parameters = parameters(velocity_dis=Normal, velocity_variance=1.0, velocity=5.0)
velocity_distribution = create_velocity_dis(parameters)
```
"""
function create_veloctiy_dis(parameters::parameters)
    distribution = parameters.pa_ph.velocity_dis
    var = parameters.velocity_variance
    velocity = parameters.velocity
    if velocity == 0 && var == 0 
        concrete_distribution = Normal(0.0,0.0) 
    else
        if distribution == Uniform
        concrete_distribution = truncated(distribution(velocity-var,velocity+var),lower = 0)
        end

        if distribution == Normal
        concrete_distribution = truncated(distribution(velocity,var),lower = 0)
        end
    end
    return concrete_distribution
end

"""
    periodic_boundary_correction(x_pos::Int64, y_pos::Int64, size::Tuple{Int64,Int64}) -> Tuple{Int64, Int64}

Corrects the positions `(x_pos, y_pos)` of an object within a 2D grid to ensure they fall within the grid's boundaries,
assuming periodic boundary conditions. This means that if an object moves beyond one edge of the grid, it re-enters the grid
from the opposite edge, similar to a toroidal (doughnut-shaped) space.

# Arguments
- `x_pos::Int64`: The x-coordinate of the object's position.
- `y_pos::Int64`: The y-coordinate of the object's position.
- `size::Tuple{Int64,Int64}`: A tuple representing the size of the grid, where the first element is the width (x dimension) and the second element is the height (y dimension).

# Returns
- `Tuple{Int64, Int64}`: The corrected position of the object as a tuple of `(x_pos, y_pos)`.

# Examples
```julia
# Given a grid of size 10x10
grid_size = (10, 10)

# An object at position (11, 5) would wrap around to (1, 5)
periodic_boundary_correction(11, 5, grid_size) # returns (1, 5)

# An object at position (0, 10) would wrap around to (10, 10)
periodic_boundary_correction(0, 10, grid_size) # returns (10, 10)
```
"""
function periodic_boundary_correction(x_pos::Int64, y_pos::Int64, size::Tuple{Int64,Int64})
   
    x_size = size[1]  # Width of the grid
    y_size = size[2]  # Height of the grid    
    
    # Correct x_pos if it moves beyond the right edge of the grid
    if x_pos > x_size
        x_pos = x_pos - x_size
    end
    # Correct x_pos if it moves beyond the left edge of the grid
    if x_pos < 1
        x_pos = x_pos + x_size 
    end
    
    # Correct y_pos if it moves beyond the top edge of the grid
    if y_pos > y_size
        y_pos = y_pos - y_size
    end
    # Correct y_pos if it moves beyond the bottom edge of the grid
    if y_pos < 1
        y_pos = y_pos + y_size 
    end
    
    return(x_pos, y_pos)  # Return the corrected position
end
"""
- **Function Name**: `path_tracing`
- **Description**: This function simulates the movement of an agent (`tryp`) through a grid, step by step, to determine the point at which it first encounters a wall. It returns the X and Y positions on the grid where the first wall is encountered, along with the distance already walked.
- **Parameters**:
    - `grid::Matrix{Int64}`: The grid through which the agent moves, represented as a matrix of integers.
    - `tryp::agent`: The agent attempting to move through the grid. It contains properties like direction (`dir_x`, `dir_y`) and position (`x_pos`, `y_pos`).
    - `velocity`: The speed at which the agent attempts to move through the grid.
    - `parameters::parameters`: Additional parameters affecting the simulation, such as boundary conditions.
- **Returns**: The function returns a tuple containing the X and Y indices of the first wall grid point encountered and the remaining distance the agent can walk.
- **Algorithm**:
    1. Initialize the agent's direction and position.
    2. Iterate through the steps based on the agent's velocity.
    3. Update the agent's position at each step, rounding to the nearest grid point.
    4. If periodic boundary conditions are specified, adjust the position accordingly.
    5. Check if the current grid point is a wall (`grid[Y_step_i, X_step_i] < 0`). If so, move the agent back to the previous position and adjust for periodic boundary conditions if necessary.
    6. Return the final position and the distance walked before encountering the wall.
"""
function path_tracing(grid::Matrix{Int64},tryp::agent,velocity,parameters::parameters)
    dir_x = tryp.dir_x
    dir_y = tryp.dir_y
    X_step = tryp.x_pos
    Y_step = tryp.y_pos
    X_step_i = round(Int64,X_step)
    Y_step_i = round(Int64,Y_step)
    walked_dist = 0 

    for len in 1:velocity
        walked_dist = len
        #needs no normalization as dir_x and y are already normalized to 1 
        X_step += dir_x
        Y_step += dir_y
        #make integer versions of current position which are rounded to next gridpoint
        X_step_i = round(Int64,X_step)
        Y_step_i = round(Int64,Y_step)
        
        if parameters.boundary_conditions == "periodic"
            X_step_i,Y_step_i = periodic_boundary_correction(X_step_i, Y_step_i, parameters.size)
        end
        #check if gridpoint is occupied with boundary and return previous gridpoint
        if grid[Y_step_i,X_step_i] < 0
            X_step -= dir_x
            Y_step -= dir_y
            #make integer versions of current position which are rounded to next gridpoint
            X_step_i = round(Int64,X_step)
            Y_step_i = round(Int64,Y_step)
            walked_dist = len-1
            # add additional periodic check other ise move back can move outside of grid if a boundary check has performaned before
            if parameters.boundary_conditions == "periodic"
                X_step_i,Y_step_i = periodic_boundary_correction(X_step_i, Y_step_i, parameters.size)
            end
            break
        end
    end
    
    return X_step_i ,Y_step_i, walked_dist
end


"""
Returns all grid points surrounding the collision point in a radius of r 
The discretization going from the continuous circle to the discrete points on the grid is done with the stepwidth:  `steps` 
"""
function find_surounding_points(grid::Matrix{Int64},tryp::agent, Parameters::parameters)
    r = Parameters.pa_ph.radius_tanget
    steps = 8*r
    y_list =Vector{Int64}(undef,steps)
    x_list = Vector{Int64}(undef,steps)
    value_list = Vector{Int64}(undef,steps)
    #println("x = $x y = $y ")

    for i in 1:steps
        y_t = tryp.y_pos + round(Int64,sin(2*π*i/steps)*r)
        x_t = tryp.x_pos + round(Int64,cos(2*π*i/steps)*r)
        y_list[i] = y_t
        x_list[i] = x_t
        
        #add periodic boundary correction
        if Parameters.boundary_conditions == "periodic"
            x_t, y_t = periodic_boundary_correction(x_t, y_t, Parameters.size)
        end

        value_list[i] = grid[y_t,x_t]
        # binarize value list to -1 and 0 for boundary and occupied gridpositions 
        value_list[i] >= 0 ? value_list[i] = 0 : value_list[i] = -1   
    end
    return (value_list,x_list, y_list)
end


"""
Finds and returns the two most distant tangent points from a given array of tangent points.

This function is designed to process an array of tangent points, each represented as a `Vector{Int64}` containing two elements (coordinates). It calculates the pairwise distances between all points in the array using a helper function `length_vec` to determine the length of the vector formed by two points. The function then identifies and returns the pair of points that are the most distant from each other, which are considered the "right" tangent points for further processing.

# Arguments
- `tangent_points::Vector{Vector{Int64}}`: An array of points, each point is a vector of two integers representing its coordinates.

# Returns
- `Vector{Vector{Int64}}`: An array containing the two most distant tangent points from the input array.

# Example
```julia
tangent_points = [[1, 2], [3, 4], [5, 6]]
right_tangent_points = find_right_tangent_points(tangent_points)
# right_tangent_points will contain the two points that are furthest apart
```
"""
function find_right_tangent_points(tangent_points::Vector{Vector{Int64}})
    right_tangent_points = [[1, 0], [2, 0]]
    right_length = 0.0
    for i in 1:(length(tangent_points)-1)
        for j = (i+1):length(tangent_points)
            temp_vec = tangent_points[j]-tangent_points[i]
            if right_length < length_vec(temp_vec)
                right_length = length_vec(temp_vec)
                right_tangent_points = [tangent_points[j],tangent_points[i]]
            end
        end
    end
    return right_tangent_points
end

  
        
"""
Reurn the gridpoints from the list of gridpoints in circle, which are the first ones of the boundary 
by checking when there is a change in the value list of points from -1 to 0 or vice versa. The last gridpoint that has the 
value -1 will be used as the probable point on the boundary between occupied and free space. 
"""
function find_tangent_points(value_list::Vector{Int64}, x_list::Vector{Int64}, y_list::Vector{Int64})
    tangent_points = Vector{Vector{Int64}}(undef,0)
    # compare neigbouring array elements
    for i in 1:length(value_list)
        # account for array overflow
        j = i+1
        j > length(value_list) && (j = 1)
        # asing border point to index with a negative value
        if value_list[i] != value_list[j]
            if value_list[i] >= 0
                push!(tangent_points, [x_list[j],y_list[j]])
            else
                push!(tangent_points, [x_list[i],y_list[i]])
            end
        end
    end
    tangent_points = unique(tangent_points)
    if length(tangent_points) < 2
        #println("less  than two tangent points found :$(tangent_points), no tangent defineable")
    end

    if length(tangent_points) > 2

        tangent_points = find_right_tangent_points(tangent_points)
       
        #println(length(tangent_points))
        #throw("More than two tangent point found, no unambiguous tangent defined")
    end
    return tangent_points
end


function rotate_vec(vec::Vector{Float64},Θ::Float64) 
    rot_mat = [cos(Θ) -sin(Θ); sin(Θ) cos(Θ)]
    rot_vec = rot_mat*vec
    return rot_vec
end


"""
Calculates outgoing path of incoming particle, using raytracing formular and generic conventions: 
#![reflection_calculus](https://slidetodoc.com/presentation_image_h/ea91543057889988b7fa9402973b147a/image-11.jpg)
"""
function calculate_outgoing_vec(tryp::agent, normal_vec::Vector{Float64})
    l = -1* [tryp.dir_x,tryp.dir_y]
    n = normal_vec
    r = 2*(l ⋅ n)*n -l
    return r
end


"""
Return the normalized reflected vector from the tangent to the next boundary of a given grid point. Function includes a fallback 
option, if no unambiguous tangent is defined by two points, the original orientation is just reversed.   
"""
function get_reflected_vec(grid::Matrix{Int64}, tryp::agent,Parameters::parameters)
    
    a,b,c =find_surounding_points(grid::Matrix{Int64},tryp::agent,Parameters::parameters)
    tangent_points = find_tangent_points(a,b,c)
    # if tangent point calculation has failed 

    # failback to just reverseing original orientation 
    if length(tangent_points) != 2
        reflected_vec = -1*[tryp.dir_x, tryp.dir_y]
        #global x += 1
        #println("reflection_vec_calc_failed_orientation_just reversed ", x)
    else
        tangent_vec = tangent_points[1]-tangent_points[2]
        #cast to float otherwise normalize can cause errors
        tangent_vec = normalize(Float64.(tangent_vec))
        normal_vec = rotate_vec(tangent_vec, Float64(π)/2)
        #normal_vec = normalize(normal_vec)
        reflected_vec = calculate_outgoing_vec(tryp,normal_vec)
    end
    return reflected_vec
end



""" 
Calculates reflection of agent on wall, calls itself up to 2 times recursivley if multiple reflections occur

#![reflection_calculus](https://slidetodoc.com/presentation_image_h/ea91543057889988b7fa9402973b147a/image-11.jpg)
"""
function calc_path(grid::Matrix{Int64}, fourier_grid::Matrix{Int},tryp::agent, velocity, Parameters::parameters, counter = 1)
    X_col, Y_col, Walked_dist = path_tracing(grid,tryp,velocity,Parameters)
    leftover_dist = velocity-Walked_dist
    # create new tryp which position is directly on the border or at the end of the traced path
    # this agent can either be returned or used as a helping datastructure in case of reflections
    tryp_temp = agent(X_col, Y_col, tryp.dir_x, tryp.dir_y)
    if leftover_dist == 0
        return tryp_temp

    else
        reflected_vec = get_reflected_vec(grid,tryp_temp,Parameters)
        Dir_x_new = reflected_vec[1]
        Dir_y_new = reflected_vec[2]
        
        if Parameters.boundary_moveable == true
            remove_boundary!(grid,fourier_grid,tryp_temp,Parameters)
        end
        counter += 1
        # condition that prevent infinite loop, so that only 3 consecutive reflections are possible 
        if counter >4
            #println(counter)
            return agent(X_col,Y_col,Dir_x_new,Dir_y_new)
        else
            return calc_path(grid,fourier_grid, agent(X_col,Y_col,Dir_x_new,Dir_y_new), leftover_dist, Parameters,counter) 
        end
    end
end

"""
    remove_boundary!(grid::Matrix{Int64}, tryp::agent, Parameters::parameters)

Modify `grid` in-place by incrementing the value of cells within a specified radius around the position of `tryp` if those cells are part of the unmoveable space meaning that they, 
 have negative values. This operation remove/weakens the unmoveable space boundary around the agent's position by increcments of 1.

# Arguments
- `grid::Matrix{Int64}`: A 2D grid (matrix) where each cell's value represents the state of that cell. Negative values indicate unmoveable space.
- `tryp::agent`: An object representing an agent, containing fields `x_pos` and `y_pos` which denote the agent's current position on the grid.
- `Parameters::parameters`: An object containing simulation parameters, including `radius_collision` which defines the radius around the agent within which the function operates,
 and `size` for boundary corrections.

# Operation
Iterates over a square area centered on the agent's position, with side lengths of `2*radius + 1`.
 For each cell within this area, if the cell is within a circular region defined by `radius_collision`, the function checks if the cell's value is negative.
  If so, the value is incremented by 1, after applying periodic boundary corrections to account for the grid's edges.
"""

function remove_boundary!(grid::Matrix{Int64},fourier_grid::Matrix{Int},tryp::agent,Parameters::parameters)
    x_pos = tryp.x_pos
    y_pos = tryp.y_pos
    radius = Parameters.pa_ph.radius_collision
    diamond_fac = 0.10
    radius_int = radius + round(Int, radius*diamond_fac+1)
    radius = radius +  radius*diamond_fac
    for y in -radius_int:1:radius_int
        for x in -radius_int:1:radius_int
            if sqrt(x^2+y^2) <= Parameters.pa_ph.radius_collision || (abs(x) + abs(y)) < radius
        
            #if (abs(x)^(1/2) + abs(y)^(1/2))^2 < radius*2 
            #if sqrt(x^2+y^2) < radius  && (abs(x^4) + abs(y^4))^(1/4) < radius* .8
            #if (abs(x) + abs(y)) < radius
                y_i = y + y_pos
                x_i = x + x_pos

                x_i_c, y_i_c = periodic_boundary_correction(x_i,y_i, Parameters.size)
                
                
                if grid[  y_i_c,x_i_c ] < 0
                    #grid[ y_i_c,x_i_c, ] += (round(Int,(abs(fourier_grid[ y_i_c,x_i_c, ])^.6))+1)
                    grid[ y_i_c,x_i_c, ] += 1
                    if grid[ y_i_c,x_i_c, ] > 0
                        grid[ y_i_c,x_i_c, ] = 0
                    end
                end
                """
                if grid[  y_i_c,x_i_c ] < 0
                    k = sum((grid[y_ic-1:y_i+1,x_i_c-1:x_i_c+1] .< 0) .* [0 1 0 ; 1 0 1; 0 1 0 ])
                    if k > 0
                        grid[ y_i_c,x_i_c, ] += k^2
                        if grid[ y_i_c,x_i_c, ] > 0
                            grid[ y_i_c,x_i_c, ] = 0
                        end
                    else
                    grid[ y_i_c,x_i_c, ] += 1
                    end
                end
                """
            end
        end
    end
end


"""
    find_free_neighbour_points(grid::Matrix{Int64}, tryp::agent, Parameters::parameters) -> Tuple{Bool, Int64, Int64}

Finds a free neighboring point around a given agent's position in a grid, considering specified boundary conditions.

# Arguments
- `grid::Matrix{Int64}`: A 2D grid representing the environment, where 0 indicates a free space and other values indicate occupied spaces.
- `tryp::agent`: An agent object with fields `x_pos` and `y_pos` indicating its current position in the grid.
- `Parameters::parameters`: A parameters object containing simulation parameters, including `boundary_conditions` which can be "periodic" and `size` indicating the size of the grid.

# Returns
- `sucess::Bool`: A boolean indicating whether a free neighboring point was found.
- `X_pos_new::Int64`: The x-coordinate of the found free neighboring point. Returns the agent's original x-coordinate if no free point is found.
- `Y_pos_new::Int64`: The y-coordinate of the found free neighboring point. Returns the agent's original y-coordinate if no free point is found.

# Description
This function attempts to find a free (unoccupied) neighboring point around the agent's current position. It randomly shuffles the directions to check for free spaces to ensure unbiased movement direction. If the `boundary_conditions` are set to "periodic", it will apply periodic boundary corrections to wrap around the grid edges. The function returns whether a free space was found and the coordinates of that space.
"""
function find_free_neighbour_points(grid::Matrix{Int64},tryp::agent,Parameters::parameters)
    X_pos = tryp.x_pos
    Y_pos = tryp.y_pos
    X_pos_new = tryp.x_pos
    Y_pos_new = tryp.y_pos
    perm_vec = [[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,-1],[-1,1],[1,-1]]
    sucess = false
    r.shuffle!(perm_vec)
    x_perm =[i[1] for i in perm_vec]
    y_perm =[i[2] for i in perm_vec]
    for j in 1:length(x_perm)
        Y_pos_new = Y_pos+y_perm[j]
        X_pos_new = X_pos+x_perm[j]
        if Parameters.boundary_conditions == "periodic"
            X_pos_new, Y_pos_new = periodic_boundary_correction(X_pos_new,Y_pos_new,Parameters.size)
        end
        if grid[Y_pos_new, X_pos_new] == 0
            sucess = true
            break 
        end
    end
    return sucess, X_pos_new, Y_pos_new
end


function conv(img::AbstractArray, kernel::AbstractArray)
    # Calculate the size of the convolution
    csize = size(img) .+ size(kernel) .- 1

    # Initialize padded versions of the image and kernel
    padimg = zeros(ComplexF32, csize)
    padkernel = zeros(ComplexF32, csize)

    # Copy the image and kernel into the padded versions
    padimg[1:size(img,1), 1:size(img,2)] .= img
    padkernel[1:size(kernel,1), 1:size(kernel,2)] .= kernel

    # Perform the Fast Fourier Transform on the padded image and kernel
    fft!(padimg)
    fft!(padkernel)

    # Multiply the transformed image and kernel
    padimg .*= padkernel

    # Perform the inverse Fast Fourier Transform on the result
    ifft!(padimg)

    # Extract the real part of the result
    output = real.(padimg)

    # Calculate the offset for the kernel
    off = div.(size(kernel), 2)

    # Extract the convolved image from the output
    h, w = size(img)
    return output[1+off[1]:h+off[1], 1+off[2]:w+off[2]]
end

"""
    update_position!(grid::Matrix{Int64}, full_agent_list, Parameters::parameters, timestep=1)

Updates the positions of agents on a grid based on their velocities and directions.

## Arguments
- `grid::Matrix{Int64}`: A 2D grid representing the environment where each cell can be empty (0), occupied by an agent (positive non-zero) or part of unmoveable space (negative non-zero).
- `full_agent_list`: A collection of agents. Can be a `CircularBuffer{Vector{agent}}` for interactive/infinite simulations or a simple `Vector{agent}` for static/finite simulations.
- `Parameters::parameters`: A struct containing simulation parameters, including velocity distribution, boundary conditions, and whether path tracing is enabled.
- `timestep=1`: The current timestep of the simulation. Defaults to 1.

## Behavior
- Determines the velocity distribution for agents based on the provided parameters.
- Draws a random velocity for each agent from the distribution. 
- Decides on the agent list to use based on the type of simulation (interactive/infinite vs. static/finite).
- Iterates through each agent, calculating their new position based on their velocity and direction.
- If `path_tracing` is disabled, agents move directly to their new position if it's empty, respecting boundary conditions.
- If `path_tracing` is enabled, calculates a path for the agent and moves them along this path and interacts with the non movable space,
     also respecting boundary conditions and handling sliding boundary conditions if enabled.
- Updates the grid to reflect the new positions of the agents.
- Updates the non movable space if the boundary is moveable.

## Note
- This function modifies the `grid` and `full_agent_list` in place.
- Agents are assumed to be immutable; thus, to update an agent's position, a new agent with the updated position is created and replaced in the list.
"""
function update_position!(Grids::Tuple{Matrix{Int64}, Matrix{Float64}}, full_agent_list, Parameters::parameters,timestep =1)
    grid = Grids[1]
    """
    fourier_grid = Grids[1] .< 0
    laplac_kernel = [0 1 0; 1 -4 1; 0 1 0]
    fourier_grid = round.(Int, conv(fourier_grid, laplac_kernel))
    """
    fourier_grid = ones(Int, 2,2)
    velocity_distribution = create_veloctiy_dis(Parameters)
    boundary_conditions = Parameters.boundary_conditions
    #using next timestep instead of previous one as the direction and previous position 
    # have already been written there in the update direction function 

    # decision between interactive/infinte and static/finite simulation
    if typeof(full_agent_list) == CircularBuffer{Vector{agent}}
        agent_list = full_agent_list[end]
    else
        agent_list = full_agent_list[timestep+1]
    end
    # agent in agent_list does not work as agent itself cannot be modifed as it is immutable
    for i = 1:length(agent_list)
        X_pos = agent_list[i].x_pos
        Y_pos = agent_list[i].y_pos
        Dir_x = agent_list[i].dir_x
        Dir_y = agent_list[i].dir_y
        velocity = rand(velocity_distribution)

        ## without pathtracing much more simple update sheme can be used 
        if Parameters.path_tracing == false 
            #periodic boundary_conditions
            X_pos_new = X_pos + round(Int64,Dir_x*velocity)
            Y_pos_new = Y_pos + round(Int64,Dir_y*velocity)

            if boundary_conditions == "periodic"
                X_pos_new, Y_pos_new =  periodic_boundary_correction(X_pos_new, Y_pos_new, Parameters.size)
            end
            # if target gird position is empty move, else do nothing 
            if grid[Y_pos_new, X_pos_new] == 0 
                grid[Y_pos_new, X_pos_new] = grid[Y_pos, X_pos]
                grid[Y_pos, X_pos] = 0
                
                agent_list[i] = agent(X_pos_new,Y_pos_new,Dir_x,Dir_y)
            end
        elseif Parameters.path_tracing == true 
            # round velocity to next int otherwise strange stuff happens
            velocity = round(Int64, velocity)
            new_agent = calc_path(grid, fourier_grid,agent_list[i], velocity, Parameters)
            X_pos_new = new_agent.x_pos
            Y_pos_new = new_agent.y_pos

            if grid[Y_pos_new, X_pos_new] == 0
                grid[Y_pos_new, X_pos_new] = grid[Y_pos, X_pos]
                grid[Y_pos, X_pos] = 0
                agent_list[i] = new_agent
            
            elseif Parameters.sliding_movement == true
                sucess, X_pos_new, Y_pos_new = find_free_neighbour_points(grid,new_agent,Parameters)
                if sucess
                    grid[Y_pos_new, X_pos_new] = grid[Y_pos, X_pos]
                    grid[Y_pos, X_pos] = 0
                    new_agent = agent(X_pos_new,Y_pos_new,new_agent.dir_x,new_agent.dir_y)
                    agent_list[i] = new_agent
                end
            end
        end
    end
end

