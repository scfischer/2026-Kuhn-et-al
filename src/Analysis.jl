"""
    create_data_path_vector(data_path) -> Vector{String}

Creates a vector of full data paths by reading the contents of the specified directory or directories, excluding certain files.

# Arguments
- `data_path`: A string representing a single directory path or a vector of directory paths.

# Returns
- `Vector{String}`: A vector of full paths to the data files in the specified directory or directories, excluding "settings.csv", "settings.jls", and "analysis".

# Example Usage
```julia
data_path = "C:\\data"
full_data_paths = create_data_path_vector(data_path)

data_paths = ["C:\\data1", "C:\\data2"]
full_data_paths = create_data_path_vector(data_paths)
```
"""
function create_data_path_vector(data_path)
    typeof(data_path) == String ? path = [data_path] : path = data_path
    data_paths = setdiff!.(readdir.(path), fill(["settings.csv","settings.jls","parameters.jls","parameters.txt","analysis"],length(path)))
    full_data_path = [path[x] .*"\\".*data_paths[x] for x in 1:length(path)]
    return full_data_path
end


"""
    create_para_struct_vec(data_used) -> Vector{Vector{parameters}}

Creates a vector of parameter structures from the given data. The function handles both 1-dimensional and 2-dimensional data inputs.

# Arguments
- `data_used`: The input data from which parameter structures are to be created. It can be a 1-dimensional or 2-dimensional array.

# Returns
- `Vector{Vector{parameters}}`: A vector of vectors, where each inner vector contains parameter structures created from the input data.

# Example Usage
```julia
data_used = DataFrame(...)  # Assume this is a DataFrame with the necessary data
para_vec = create_para_struct_vec(data_used)
```
"""
function create_para_struct_vec(data_used)
    para_vec = Vector{Vector{parameters}}(undef,0)
    if length(size(data_used)) == 1
        para_dict = make_para_dict_list(Dict(names(data_used) .=> values(data_used)))
        push!(para_vec, create_para_struct.(para_dict))
    else
        for x in 1:size(data_used)[1]
            para_dict = make_para_dict_list(Dict(names(data_used[x,:]) .=> values(data_used[x,:])))
            push!(para_vec, create_para_struct.(para_dict))
        end
    end
    return para_vec
end


"""
    create_sweep_para_vec(data_used::DataFrame) -> Array{Array{Float64}}

Generates a vector of sweep parameter vectors based on the input data.

# Arguments
- `data_used::DataFrame`: An DataFrame where each row contains the sweep parameter information, including the sweep parameter name, start value, increment, and number of steps.

# Returns
- `Array{Array{Float64}}`: A vector of vectors, where each inner vector contains the sweep parameter values for a specific row in the input data.

# Example
```julia
data_used = Dataframe(
    (sweep_parameter="param1", param1=0.0, parameter_increment=0.1, parameter_steps=5,...),
    (sweep_parameter="param2", param2=1.0, parameter_increment=0.2, parameter_steps=3,...)
)
sweep_para_vec = create_sweep_para_vec(data_used)
# sweep_para_vec is now [[0.0, 0.1, 0.2, 0.3, 0.4], [1.0, 1.2, 1.4]]
```
"""
function create_sweep_para_vec(data_used)
    sweep_para_vec = []
    for i in 1:size(data_used)[1]
        sweep_para = data_used[i,:].sweep_parameter
        start_value = data_used[i,:][sweep_para]
        sweep_vec = [start_value + data_used[i,:].parameter_increment * x for x in 0:(data_used[i,:].parameter_steps -1)  ]
        push!(sweep_para_vec, sweep_vec)
    end
    return sweep_para_vec
end


"""
    ana_order_para_vicsek(data, para::parameters; time_steady=0.7) -> Float64

Calculates the average order parameter for the Vicsek model over the steady-state period of the simulation.

# Arguments
- `data`: A collection of data points, where `data[1]` contains the time series data for which the order parameter is to be calculated.
- `para::parameters`: A struct containing the parameters required by the `order_parameter` function.
- `time_steady::Float64=0.7` (optional): The fraction of the total time steps to be considered as the steady-state period. Default is 0.7.

# Returns
- `Float64`: The average order parameter over the steady-state period.

# Example Usage
```julia
data = [...]  # Assume this is the data collected from the simulation
para = parameters(...)  # Assume this is the parameters struct
average_order = ana_order_para_vicsek(data, para, time_steady=0.8)
```
"""
function ana_order_para_vicsek(data,para::parameters; time_steady = 0.7)
    time_steps = length(data[1])
    v_a = 0
    for i in round(Int,time_steps*time_steady):time_steps
           v_a += order_parameter(data[1][i],para)
        end
    return v_a/(time_steps-round(Int,time_steps*time_steady))
end
    

function order_parameter(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters)
    grid = grids[1]
    x_dir_sum = sum([i.dir_x for i in agent_list])
    y_dir_sum = sum([i.dir_y for i in agent_list])
    number_agents = length(agent_list)
    v_a = sqrt(x_dir_sum^2+y_dir_sum^2)*(1/number_agents)
    return v_a
end 

"""
    ana_order_para_circle(data, para::parameters; time_steady::Float64=0.7) -> Float64

Calculates the average order parameter over a specified steady-state period.

# Arguments
- `data::save_model_data`: An tuple contain an array as first element where each array element represents data points over time.
- `para::parameters`: A structure containing parameters required for the `order_parameter_circle` function.
- `time_steady::Float64`: A fraction of the total time steps to consider as the steady-state period (default is 0.7).

# Returns
- `Float64`: The average order parameter over the steady-state period.

# Example
```julia
data = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]]
para = parameters(...)
result = ana_order_para_circle(data, para)
# result is the average order parameter over the steady-state period
```
"""
function ana_order_para_circle(data,para::parameters; time_steady = 0.7)
    time_steps = length(data[1])
    v_a = 0
    for i in round(Int,time_steps*time_steady):time_steps
           v_a += order_parameter_circle(data[1][i],para)
        end
    return v_a/(time_steps-round(Int,time_steps*time_steady))
end


function order_parameter_circle(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters)
    grid = grids[1]
    number_agents = length(agent_list)
    center_x = round(Int,para.size[2]/2)
    center_y = round(Int,para.size[1]/2)
    x_pos = [i.x_pos for i in agent_list]
    y_pos = [i.y_pos for i in agent_list]

    x_lot = x_pos .- center_x
    y_lot = y_pos .- center_y
    lot_vecs = [[x_lot[i],y_lot[i]] for i in 1:number_agents]
    lot_array = normalize.(lot_vecs)
    # deleteat all lot vectors that are not definable, happens when agent is perfectly in center of system
    #index_of_false_agents = (lot_array .=== fill([NaN,NaN],length(lot_array)))
    #deleteat!(lot_array, index_of_false_agents)
    rot_mat = [0 -1; 1 0]
    rot_mat_list = fill(rot_mat, number_agents)
    tan_array = rot_mat_list.*lot_array
    index_of_false_agents = find_nan(tan_array)
    deleteat!(tan_array, index_of_false_agents)
    dir_array = [[i.dir_x,i.dir_y] for i in agent_list]
    deleteat!(dir_array, index_of_false_agents)
    v_c = sqrt((sum(dot.(dir_array,tan_array)))^2)*1/number_agents

    return v_c
end 


function order_parameter_perpendic(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters)
    grid = grids[1]
    number_agents = length(agent_list)
    center_x = round(Int,para.size[2]/2)
    center_y = round(Int,para.size[1]/2)
    x_pos = [i.x_pos for i in agent_list]
    y_pos = [i.y_pos for i in agent_list]

    x_lot = x_pos .- center_x
    y_lot = y_pos .- center_y
    lot_vecs = [[x_lot[i],y_lot[i]] for i in 1:number_agents]
    lot_array = normalize.(lot_vecs)
    # deleteat all lot vectors that are not definable, happens when agent is perfectly in center of system
    #index_of_false_agents = (lot_array .=== fill([NaN,NaN],length(lot_array)))
    #deleteat!(lot_array, index_of_false_agents)
    
    index_of_false_agents = find_nan(lot_array)
    deleteat!(lot_array, index_of_false_agents)
    dir_array = [[i.dir_x,i.dir_y] for i in agent_list]
    deleteat!(dir_array, index_of_false_agents)
    v_p = sqrt((sum(dot.(dir_array,lot_array)))^2)*1/number_agents

    return v_p
end 

"""
    radial_distribution(agent_list::Vector{agent}, size::Tuple{Int, Int})

Calculates the radial distribution of agents from the center of the grid.

# Arguments
-  agent_list::Vector{agent} : A vector of agent objects.
- `size::Tuple{Int, Int}`: A tuple containing the dimensions of the grid.

# Description
This function calculates the radial distribution of agents from the center of the grid. It performs the following steps:
1. Calculates the center of the grid.
2. Computes the relative positions of the agents with respect to the center.
3. Calculates the Euclidean distance of each agent from the center.
4. Bins the distances and normalizes by 1/(2r+1) to get the radial distribution.

# Example
```julia
x_pos = [1, 3, 5, 7]
y_pos = [2, 4, 6, 8]
size = (10, 10)

radial_distribution(x_pos, y_pos, size)
```
# Returns
- `Vector{Int}`: The radial distribution 
- `Vector{Float}`:   The distance vector, useable for histogram plotting    
"""
function radial_distribution(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters; dist_vec_return = false)
    grid = grids[1]
    size = para.size
    x_pos = [i.x_pos for i in agent_list]
    y_pos = [i.y_pos for i in agent_list]

    center_x = round(Int,size[1] ÷ 2)
    center_y = round(Int,size[2] ÷ 2)

    number_agents = length(x_pos)

    x_lot = x_pos .- center_x
    y_lot = y_pos .- center_y
    lot_vecs = [[x_lot[i],y_lot[i]] for i in 1:number_agents ]
    
    dist_vec  = norm.(lot_vecs)
    
    max_dist = floor(Int, maximum(dist_vec))+1
    dist_vec_equal = zeros(Int, max_dist)
    for i in 1:length(dist_vec)
        dist_vec_equal[floor(Int, dist_vec[i])+1] += 1
    end
    # normalize by 1/(2r+1) as the are of the circle rings increase with this factor
    dist_vec_equal = dist_vec_equal .* (1 ./ (collect(1:max_dist)*2 .-1))

    if dist_vec_return 
        return dist_vec_equal, dist_vec
    else 
        return dist_vec_equal
    end
end


function b_w_im(grid)
    grid = replace(u -> (u<0) ? 0 : 1, grid)
    return grid
end

function b_w_oc(grid)
    grid = replace(u -> (u<1) ? 0 : 1, grid)
    return grid
end

"""
    angular_metric(agent_list::Vector{agent}, grid::AbstractArray, para::parameters; steps::Int = 360)

Calculates the angular distribution of non-zero pixels in the grid relative to the center.

# Arguments
- `agent_list::Vector{agent}`: A vector containing the agents.
- `grid::AbstractArray`: The grid representing the simulation space.
- `para::parameters`: A struct containing various parameters, including the size of the grid.
- `steps::Int`: The number of angular steps to divide the circle into (default is 360).

# Description
This function calculates the angular distribution of non-zero pixels in the grid relative to the center. 
    It performs the following steps:
1. Initializes a vector for the angular metric.
2. Calculates the center of the grid.
3. Determines the size of each angular step in radians.
4. Iterates over each pixel in the grid and calculates the vector from the center to the pixel.
5. Calculates the angle of the vector and determines the corresponding sector.
6. Increments the count for the sector.
7. Adjusts the first and last sectors to account for discretization artifacts by taking their mean.
8. Returns the angular metric.
"""
function angular_metric(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters; steps::Int = 360)
    # Initialize the vector for the angular metric
    grid = grids[1] 
    angle_vec = zeros(Int64, steps)
    center = [para.size[1] ÷ 2, para.size[2] ÷ 2] 
    # Calculate the size of each step in radians
    stepsize = 2π/steps 
    start_grid = create_grids(para, create_scal_grid = false)

    grid_c = b_w_im(grid)-b_w_im(start_grid)

    # For each pixel in the image
    for y in 1:size(grid_c)[1], x in 1: size(grid_c)[2]
        # If the pixel is not zero
        if grid_c[y,x] != 0
            # Calculate the vector from the center to the pixel
            vec = [y,x] .- center

            # Calculate the angle of the vector
            ϕ = atan(vec...) + π

            # Determine the sector to which the angle belongs
            step_index = minimum([round(Int, ϕ / stepsize), steps - 1])

            # Increment the count for the sector
            angle_vec[step_index + 1] += 1
        end
    end

    # The first and last sectors may be inaccurate due to discretization artifacts, so take the mean of both
    mean_angle = round(Int, mean([angle_vec[1], angle_vec[end]]))
    angle_vec[1], angle_vec[end] = mean_angle, mean_angle

    # Return the angular metric
    return angle_vec
end



function pair_cor_metric(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters; max_samples::Int = 1000000, steps::Int = 360)
    # Initialize the vector for the pair correlation metric
    angle_vec = zeros(Int64, steps)
    grid = grids[1] 
    # Calculate the size of each step in radians
    stepsize = π/steps
    center = [para.size[1] ÷ 2, para.size[2] ÷ 2]

    start_grid = create_grids(para, create_scal_grid = false)

    grid_c = b_w_im(grid)-b_w_im(start_grid)
    # Find the indices of the non-zero pixels
    indizies = Tuple.(findall(x-> x != 0 ,grid_c))

    if length(indizies) < 30
        return angle_vec
    else 
    # For dynamic sampling of the pair correlation metric
    if binomial(length(indizies),2) < max_samples
        samples = binomial(length(indizies),2) 
    else
        samples = max_samples
    end
        for i in 1:samples
            # Randomly select two non-zero pixels
            v1, v2  = sample(indizies, 2, replace = true)
            
            # If the reference point is one of the pixels, select two new pixels
            while Tuple(center) in [v1,v2]
                v1, v2  = sample(indizies, 2, replace = true)
            end

            # Calculate the vectors from the reference point to the pixels
            v1_v = v1 .- center 
            v2_v = v2 .- center

            # Calculate the cosine of the angle between the vectors
            x = dot(v1_v,v2_v)/(norm(v1_v)*norm(v2_v))

            # Calculate the angle between the vectors
            angle = acos(x >1.0 ? 1.0 : x < -1.0 ? -1.0 : x)

            # Determine the sector to which the angle belongs
            step_index = trunc(Int64,angle/stepsize)

            # If the step index is 360 due to rounding errors, set it to 359
            step_index_final = (step_index == 360 ? 359 : step_index)

            # Increment the count for the sector
            angle_vec[step_index_final+1] += 1
        end

        # Return the pair correlation metric
        return angle_vec/mean(angle_vec)   
    end
end


function adsorbed_material(agent_list::Vector, grids::Tuple{AbstractArray,AbstractArray}, para::parameters)

    T = grids[2]

    return sum(T)
end

    
"""
    area_gain(agent_list::Vector{agent}, grids::Tuple{AbstractArray,AbstractArray}, para::parameters)

Calculate the realtive gain in area walkable by agents in the grid.

# Arguments
- `agent_list::Vector{agent}`: A list of agents.
- `grids::Tuple{AbstractArray, AbstractArray}`: A tuple containing the grids.
- `para::parameters`: The parameters for the simulation.

# Returns
- `area_gain::Float64`: The ratio of the current area covered by agents to the original area.

# Description
This function calculates the realtive gain in area covered by agents in the grid. 
    It compares the current area walkable by agents to the original area walkable by agents at the 
    start of the simulation. The original area is calculated using a newly created grid with the 
    same parameters, and the current area is calculated using the provided grid.

# Example
```julia
agent_list = [...]  # List of agents
grids = (grid1, grid2)  # Tuple of grids
para = parameters(...)  # Simulation parameters
gain = area_gain(agent_list, grids, para)
```
"""
function area_gain(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters)
    grid = grids[1]
    original_area = sum(b_w_im(create_grids(para, create_scal_grid = false)))
    current_area = sum(b_w_im(grid))
    return current_area/original_area -1
end


function fourier_ik(metric, Og_size)
    fourier = fft(metric)
    fourier = abs.(fourier[2:360÷2])
    fourier = fourier ./ Og_size
    return sum(fourier)
end

function roughness(metric, OG_size::Int)
    return std(metric/OG_size)
 end

function cell_gain(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters)
    og_agent_number = para.pa_ph.agent_number
    
    current_agents = length(agent_list)
    return current_agents/og_agent_number -1
end

"""
    radial_density_c(agent_list::Vector{agent}, grids::Tuple{AbstractArray,AbstractArray}, para::parameters)

Calculate the radial density of agents in a grid.

# Arguments
- `agent_list::Vector{agent}`: A list of agents.
- `grids::Tuple{AbstractArray, AbstractArray}`: A tuple containing the grids.
- `para::parameters`: The parameters for the simulation.

# Returns
- `dist_vec_equal_oc::Vector{Float64}`: A vector representing the normalized radial density of agents.

# Description
This function calculates the radial density of agents in a grid by comparing all possible available
     grid points where agents can be with the points where they actually are.
      It sorts these points by their distance from the center of the grid and normalizes the 
      counts of these distances to provide a radial density distribution.

# Example
```julia
agent_list = [...]  # List of agents
grids = (grid1, grid2)  # Tuple of grids
para = parameters(...)  # Simulation parameters
density = radial_density_c(agent_list, grids, para)
```
"""
function radial_density_c(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters)
    
    data_oc = b_w_oc(grids[1])
    data_em = b_w_im(grids[1]) 
    
    center = CartesianIndex(round.(Int, size(grids[1])./2))

    lot_vec_oc = Tuple.(findall(x -> x == 1, data_oc) .-center)
    dist_vec_oc = norm.(lot_vec_oc)

    lot_vec_em = Tuple.(findall(x -> x == 1, data_em ).-center)
    dist_vec_em = norm.(lot_vec_em)

    max_dist = floor(Int, maximum(dist_vec_em))+1

    dist_vec_equal_oc = zeros(Int, max_dist)
    dist_vec_equal_em = zeros(Int, max_dist)
    for i in 1:length(dist_vec_em)
        if i<length(dist_vec_oc)
            dist_vec_equal_oc[floor(Int, dist_vec_oc[i])+1] += 1
        end
        dist_vec_equal_em[floor(Int, dist_vec_em[i])+1] += 1
    end

    return dist_vec_equal_oc./dist_vec_equal_em
end


"""
    moving_avg(X::Vector, numofele::Int=length(X)÷2)

Calculates the moving average of a vector.

# Arguments
- `X::Vector`: The input vector for which the moving average is to be calculated.
- `numofele::Int`: The number of elements to consider for the moving average (default is half the length of `X`).

# Description
This function calculates the moving average of the input vector `X`. It performs the following steps:
1. Determines the number of elements to look behind and ahead for the moving average.
2. Initializes an output vector `Y` of the same length as `X`.
3. Iterates over each element in `X` and calculates the moving average for a range of elements around it.
4. Stores the calculated moving average in the output vector `Y`.

# Example
```julia
X = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
numofele = 3

moving_avg(X, numofele)
```
# Returns
- `Y::Vector`: The vector containing the moving average values.
"""
function moving_avg(X::Vector,numofele::Int= length(X)÷2)
    # Calculate the number of elements to look behind and ahead for the moving average
    BackDelta = div(numofele,2) 
    ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
    
    len = length(X)
    
    # Initialize the output vector
    Y = zeros(Float64,len)
    
    # Calculate the moving average for each element in the input vector
    for n = 1:len
        # Determine the range of elements to include in the moving average
        lo = max(1,n - BackDelta)
        hi = min(len,n + ForwardDelta)
        
        # Calculate the moving average and store it in the output vector
        Y[n] = mean(X[lo:hi])
    end
    return Y
end



function find_nan(Lot_vecs)
    index_vec = Int[]
    for (i,x) in enumerate(Lot_vecs)
        if true in isnan.(x)
            push!(index_vec, i)
        end
    end
    return index_vec
end


"""
    nearby_neighbours(Agent_list::Vector{agent}, Grid::AbstractArray, Parameters::parameters)

Legacy function to calculate the number of nearby neighbours for each agent within a 
    specified interaction range as an heatmap to then calculate the 
    average nearby topology of an agent.

# Arguments
- `Agent_list::Vector{agent}`: A vector containing the agents.
- `Grid::AbstractArray`: The grid representing the simulation space.
- `Parameters::parameters`: A struct containing various parameters, including the interaction radius.

# Description
This function calculates the number of nearby neighbours for each agent within a specified interaction range. It performs the following steps:
1. Initializes a sum grid to accumulate the number of neighbours.
2. Iterates over each agent and extracts a sub-grid around the agent's position within the interaction range.
3. Replaces positive values in the sub-grid with 1 and negative values with 0.
4. Sets the agent's own position in the sub-grid to 0.
5. Adds the sub-grid to the sum grid.
6. Returns the rotated sum grid.

# Example
```julia
Agent_list = [agent(x_pos=5, y_pos=5), agent(x_pos=10, y_pos=10)]
Grid = rand(-1:1, 20, 20)
Parameters = parameters(interaction_radius=2)

nearby_neighbours(Agent_list, Grid, Parameters)
Returns
sum_grid: A grid representing the number of nearby neighbours for each agent. 
"""
function nearby_neighbours(Agent_list,Grid,Parameters)
    interaction_range = Parameters.interaction_radius
    sum_grid = zeros(Int64,2*interaction_range+1,2*interaction_range+1)
    x_size = size(Grid)[2]
    y_size = size(Grid)[1]
    for (i,agent) in enumerate(Agent_list)
        x_start = agent.x_pos-interaction_range
        x_end = agent.x_pos+interaction_range
        y_start = agent.y_pos-interaction_range
        y_end = agent.y_pos+interaction_range
        #create slice of grid in in interaction range
        if x_start >= 1 && x_end < x_size && y_start >= 1 && y_end < y_size
            sub_grid = Grid[y_start:y_end,x_start:x_end]
            replace!(x -> (x>0) ? 1 : x, sub_grid)
            replace!(x -> (x<0) ? 0 : x, sub_grid)
            sub_grid[interaction_range+1,interaction_range+1] = 0
        
            sum_grid += sub_grid
        end
        
    end
    return rotr90(sum_grid)
end



function push_df!(df, settings, date, name, group = "none", path = "none")
    push!(df, (settings, date, name, group))
end



function transform_data(data; i = 0)
    if i == 0
        return [[data[1][j],(data[2][j],data[3][j]),data[4]] for j in 1:length(data[1])]
    end

    if i != 0
        return [[data[1][i],(data[2][i],data[3][i]),data[4]]]
    end
    return data
end


function radial_density(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters)
end 



function og_size(agent_list::Vector{agent},grids::Tuple{AbstractArray,AbstractArray},para::parameters)

    start_grid = create_grids(para, create_scal_grid = false)

    sum(b_w_im(start_grid))
    return sum(b_w_im(start_grid))
end