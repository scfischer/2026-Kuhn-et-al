
"""
    `create_noise_dis(parameters::parameters)`

Create a noise distribution object based on the simulation parameters.

# Arguments
- `parameters::parameters`: A struct containing the simulation parameters. This must include:
  - `noise_dis`: The type of distribution to use for noise. Expected to be either `Uniform` or `Normal`.
  - `noise_strength`: A parameter (η) that defines the strength of the noise. For a `Uniform` distribution, it defines the range `[-η, η]`. 
        For a `Normal` distribution, it defines the standard deviation σ with a mean of 0.

# Returns
- A distribution object corresponding to the specified noise distribution and strength. If `noise_strength` is 0, 
it returns a `Normal` distribution with mean 0 and standard deviation 0.

# Example
```julia
parameters = parameters(noise_dis=Uniform, noise_strength=0.5)
noise_distribution = create_noise_dis(parameters)
This will create a Uniform(-0.5, 0.5) distribution for noise.
```

# Notes
The function supports creating Uniform and Normal distributions based on the `noise_dis` parameter in `parameters`.
If `noise_strength` (η) is set to 0, the function defaults to creating a Normal(0.0, 0.0) distribution,
 effectively adding no noise.
 """
function create_noise_dis(parameters::parameters)
    distribution = parameters.pa_ph.noise_dis
    η = parameters.pa_ph.noise_strength
    if η == 0
        concrete_distribution = Normal(0.0,0.0) 
    else
        if distribution == Uniform
            concrete_distribution = distribution(-η,η)
        end

        if distribution == Normal
            concrete_distribution = distribution(0.0,η)
        end
    end
    return concrete_distribution
end


function check_interactivity(full_agent_list, parameters::parameters,timestep)
    # decision between interactive/infinte and static/finite simulation
    if typeof(full_agent_list) == CircularBuffer{Vector{agent}}
        agent_list = full_agent_list[end]
        # use the @view macro to create refernces to make code identically for both versions
        agent_list_new = deepcopy(full_agent_list[end])
        #agent_list_new = @view agent_list_new_static
    else
        agent_list = full_agent_list[timestep]

        # overwrite initialized agent_list in case of dynamic simulation with new agent_list and correct size
        if parameters.growth_rate > 0
            full_agent_list[timestep+1] = deepcopy(agent_list)
        end
        
        # create a reference to sub array in order to use same variable names as in circular version
        agent_list_new = @view full_agent_list[timestep+1][:]
    end
    return agent_list, agent_list_new
end


function add_noise_angle(new_dir::AbstractArray, parameters::parameters)
    # add random angle with 2D rotation matrix
    noise_distribution = create_noise_dis(parameters)
    η = rand(noise_distribution)*parameters.pa_ph.Δt_walker

    rot_mat = [cos(η) -sin(η); sin(η) cos(η)]
    new_dir = rot_mat*new_dir
    # error handling because in very very rare cases 
    # new dir can be no legal float, no idear why, try everythig to debug this, i failed 
    if isnan(new_dir[1]) || isnan(new_dir[2]) 
        println("Direction function __new_dir =Nan,set to 0.0")
        println(new_dir[1],new_dir[2],rot_mat,η)
        new_dir[1] = 1.0
        new_dir[2] = 0.0
    end

    return new_dir
end


"""
    update_directions!(grid::Matrix{Int64}, full_agent_list, parameters::parameters, scal_grid::Matrix{Float64}=zeros(2,2), timestep=1)

Updates the directions of agents in a simulation based on the Vicsek alignment rule, repulsion, or chemotaxis. 
The function supports both interactive (infinite) and static (finite) simulations.

# Arguments
- `grid::Matrix{Int64}`: A grid representing the environment, where each cell contains the index of an agent or 0 if the cell is empty.
- `full_agent_list`: A list or circular buffer of agent states over time.
- `parameters::parameters`: A struct containing the simulation parameters, including agent number, interaction radius, repulsion settings, and chemotaxis flag.
- `scal_grid::Matrix{Float64}=zeros(2,2)` (optional): A scalar grid used for chemotaxis calculations. Default is a dummy 2x2 matrix of zeros.
- `timestep::Int=1` (optional): The current timestep of the simulation. Default is 1.

# Returns
- `Nothing`: This function modifies the `full_agent_list` and `grid` in place.

# Example Usage
```julia
grid = zeros(Int64, 100, 100)  # Example grid
full_agent_list = CircularBuffer{Vector{agent}}(10)  # Example circular buffer of agents
parameters = parameters(agent_number=100, interaction_radius=5.0, repulsion_flag=true, repulsion_range=2.0, chemotaxis_flag=true)
scal_grid = create_grad_grid(parameters)  # Example scalar grid for chemotaxis
update_directions!(grid, full_agent_list, parameters, scal_grid, timestep=1)
```
"""

#update function with  \detla t ..... 


function update_directions!(grids::Tuple{Matrix{Int64}, Matrix{Float64}}, full_agent_list, parameters::parameters,timestep = 1)
    # decision between interactive/infinte and static/finite simulation
    agent_list, agent_list_new = check_interactivity(full_agent_list, parameters,timestep)

    acces_vector = r.shuffle([1:length(agent_list);])
    grid = grids[1]
    scal_grid = grids[2]

    if parameters.chemotaxis_flag
        grad_x, grad_y = imgradients(scal_grid,KernelFactors.sobel )
        biggest_grad = maximum(sqrt.(grad_x.^2 + grad_y.^2))
    end

    # iterate over all agents in random order to calculate new direction, 
    # based on Vicsek alingment rule,  repulsion or chemotaxis 
    Threads.@threads for i in acces_vector
    #for i in acces_vector
        X_pos = agent_list[i].x_pos
        Y_pos = agent_list[i].y_pos
        X_dir = agent_list[i].dir_x
        Y_dir = agent_list[i].dir_y



        if parameters.repulsion_flag || parameters.chemotaxis_flag == false
            all_list = find_neighbour_agents(parameters, agent_list,i, grid)
            dir_x_list = all_list[1]
            dir_y_list = all_list[2]
            dist_list = all_list[3]
        end


        # decision which interaction rule to use
        if parameters.repulsion_flag
            # find all agents in repulsion range except self, therefore > 0.5
            index_rep = findall(x->x <= parameters.repulsion_range && x > 0.1,dist_list)
            if length(index_rep) > 0
                new_dir = calc_repulsion(all_list, index_rep, agent_list[i])
            else
                if parameters.chemotaxis_flag
                    new_dir = chemotactic_alignment(grad_x,grad_y,Y_pos,X_pos,X_dir,Y_dir,biggest_grad, parameters)
                else
                    new_dir = normalize([sum(dir_x_list),sum(dir_y_list)])
                end
            end
        
        elseif parameters.repulsion_flag == false

            if parameters.chemotaxis_flag
                new_dir = chemotactic_alignment(grad_x,grad_y,Y_pos,X_pos,X_dir,Y_dir,biggest_grad, parameters)
            else
                new_dir = normalize([sum(dir_x_list),sum(dir_y_list)])
            end
        end
        new_dir = add_noise_angle(new_dir, parameters)
        

        agent_list_new[i] = agent(X_pos,Y_pos,new_dir[1],new_dir[2])
    end
    # push circular buffer forward in case of interactive simulation 
    # for static, this is not neccessary as array is already initialized with full size
    if typeof(full_agent_list) == CircularBuffer{Vector{agent}}
        push!(full_agent_list,agent_list_new)
    end

end



"""
    chose_interaction_radius(parameters::parameters) -> Float64

Determines the interaction radius to be used based on the given parameters. If the `repulsion_flag` is set,
 it compares the `repulsion_range` with the `interaction_radius` and returns the larger of the two. 
 If the `repulsion_flag` is not set, it simply returns the `interaction_radius`.

# Arguments
- `parameters::parameters`: A struct containing the simulation parameters, including:
  - `repulsion_flag::Bool`: A flag indicating whether repulsion is considered.
  - `repulsion_range::Int`: The range of repulsion.
  - `interaction_radius::Int`: The default interaction radius.

# Returns
- `Int64`: The chosen interaction radius based on the given parameters.
```
"""
function chose_interaction_radius(parameters::parameters)
    if parameters.repulsion_flag == false
        if parameters.chemotaxis_flag == false
            return parameters.interaction_radius
        else 
            return 0 # if only chemotaxis is active, no interaction radius is needed
        end
    elseif parameters.repulsion_flag
        if parameters.chemotaxis_flag
            return parameters.repulsion_range
        
        elseif parameters.repulsion_range > parameters.interaction_radius
            return parameters.repulsion_range
        else
            return parameters.interaction_radius
        end
    end
end

"""
    find_neighbour_agents(parameters, agent_list, i, interaction_radius, grid) -> Vector{Vector{Float64}}

Finds the neighboring agents within a specified interaction radius for a given agent and returns their directional vectors, distances, and relative positions.

# Arguments
- `parameters`: A struct containing the simulation parameters, including boundary conditions and grid size.
- `agent_list`: A list of agents, where each agent has attributes such as `x_pos`, `y_pos`, `dir_x`, and `dir_y`.
- `i`: The index of the agent for which neighbors are being found.
- `interaction_radius`: The radius within which to search for neighboring agents.
- `grid`: A 2D array representing the grid, where each cell contains the index of an agent or 0 if the cell is empty.

# Returns
- `Vector{Vector{Float64}}`: A vector containing five vectors:
  - `dir_x_list`: The x-components of the direction vectors of the neighboring agents.
  - `dir_y_list`: The y-components of the direction vectors of the neighboring agents.
  - `dist_list`: The distances from the given agent to each neighboring agent.
  - `pos_x_list`: The x-coordinates of the neighboring agents relative to the given agent.
  - `pos_y_list`: The y-coordinates of the neighboring agents relative to the given agent.

# Example Usage
```julia
parameters = Parameters(boundary_conditions="periodic", size=(100, 100))
agent_list = [Agent(10, 10, 1.0, 0.0), Agent(12, 12, -1.0, 0.0)]  # Example agents
grid = zeros(Int64, 100, 100)
grid[10, 10] = 1
grid[12, 12] = 2
neighbors = find_neighbour_agents(parameters, agent_list, 1, 5, grid)

end
```
This function is useful for simulations where agents interact with their neighbors within a certain radius, such as in flocking or swarming behaviors. 
It calculates the direction vectors, distances, and relative positions of neighboring agents,
 which can be used to update the agent's direction based on the interactions.
"""
function find_neighbour_agents(parameters, agent_list,i, grid)
    
    interaction_radius = chose_interaction_radius(parameters)   
    X_pos = agent_list[i].x_pos
    Y_pos = agent_list[i].y_pos
    dir_x_list = Vector{Float64}(undef, 0)
    dir_y_list = Vector{Float64}(undef, 0)
    dist_list = Vector{Float64}(undef, 0)
    pos_x_list = Vector{Int64}(undef, 0)
    pos_y_list = Vector{Int64}(undef, 0)
    for Δy = -interaction_radius:interaction_radius

        for Δx = -interaction_radius:interaction_radius
            if sqrt(Δx^2+Δy^2)> interaction_radius
                continue
            end

            y = Y_pos + Δy

            x = X_pos + Δx 

            if parameters.boundary_conditions == "periodic"
                x,y = periodic_boundary_correction(x,y,parameters.size)
            end
            
            if grid[y,x] > 0
                index = grid[y,x]
                push!(dir_x_list, agent_list[index].dir_x)
                push!(dir_y_list, agent_list[index].dir_y)
                push!(pos_x_list, Δx)
                push!(pos_y_list, Δy)
                push!(dist_list, sqrt(Δx^2+Δy^2))
            end
        end
    end
    return [dir_x_list, dir_y_list, dist_list, pos_x_list, pos_y_list]
end

"""
    calc_repulsion(all_list, index_rep, old_agent) -> Vector{Float64}

Calculates the hard repulsion vector for an agent based on the positions of other agents in its vicinity.
     The function normalizes the repulsion vectors and sums them to determine the overall direction of repulsion.

# Arguments
- `all_list`: A list containing various attributes of all agents created by the function `find_neighbour_agents`.
     The x and y positions of the agents are expected to be at indices 4 and 5, respectively.
- `index_rep`: A list of indices representing the agents that are considered for repulsion.
- `old_agent`: An agent struct containing the current direction of the agent (`dir_x` and `dir_y`).

# Returns
- `Vector{Float64}`: A normalized vector representing the direction of repulsion. If no repulsion is detected, it returns the normalized direction of the `old_agent`.

# Example Usage
```julia
all_list = [...]  # Assume this is a list containing agent attributes
index_rep = [1, 2, 3]  # Indices of agents to consider for repulsion
old_agent = Agent(1, 2, 0.5, 0.5)  # Assume this is an agent struct with direction fields
repulsion_vector = calc_repulsion(all_list, index_rep, old_agent)
```
"""
function calc_repulsion(all_list, index_rep, old_agent)

    pos_x_list = all_list[4]
    pos_y_list = all_list[5]
    rep_vec = [0.0,0.0]
    for index in index_rep
        x_rep = pos_x_list[index]  
        y_rep = pos_y_list[index]
        rep_vec_temp = [x_rep,y_rep]
        # only add repulsion if position is not [0.0, 0.0] as normalization would fail then
        if rep_vec_temp != [0.0,0.0]
            rep_vec += normalize(-[x_rep,y_rep])
        end
    end
    if rep_vec == [0.0,0.0]
        return normalize([old_agent.dir_x, old_agent.dir_y])
    else
        return normalize(rep_vec)
    end
end

"""
    cross_2d(v1::Vector{<:Real}, v2::Vector{<:Real}) -> Real

Calculates the 2D cross product (also known as the perp product) of two 2D vectors.

# Arguments
- `v1::Vector{Real}`: The first 2D vector.
- `v2::Vector{Real}`: The second 2D vector.

# Returns
- `Real`: The scalar result of the 2D cross product.

# Example
```julia
v1 = [1.0, 2.0]
v2 = [3.0, 4.0]
result = cross_2d(v1, v2)
# result is -2.0
```
"""
function cross_2d(v1::Vector{T} where T<:Real, v2::Vector{T} where T<:Real)   
    return v1[1]*v2[2] - v1[2]*v2[1]
end


"""
    angle_between_vec(v1::Vector{<:Real}, v2::Vector{<:Real}) -> Real

Calculates the signed angle between two 2D vectors `v1` and `v2` using the arctangent of the 2D cross product and the dot product.

# Arguments
- `v1::Vector{Real}`: The first 2D vector.
- `v2::Vector{Real}`: The second 2D vector.

# Returns
- `Real`: The angle between the two vectors in radians.

# Example
```julia
v1 = [1.0, 0.0]
v2 = [0.0, 1.0]
angle = angle_between_vec(v1, v2)
# angle is π/2 (1.5707963267948966 radians)
```
"""
function angle_between_vec(v1::Vector{<:Real}, v2::Vector{<:Real})
    return atan(cross_2d(v1,v2),dot(v1,v2))
end

"""
    chemotactic_alignment(grad_x, grad_y, Y_pos, X_pos, X_dir, Y_dir, biggest_grad, para::parameters)

Aligns the direction of an agent based on the chemotactic gradient.

# Arguments
- `grad_x`: The gradient of the chemical concentration in the x direction.
- `grad_y`: The gradient of the chemical concentration in the y direction.
- `Y_pos`: The y position of the agent.
- `X_pos`: The x position of the agent.
- `X_dir`: The x component of the agent's direction.
- `Y_dir`: The y component of the agent's direction.
- `biggest_grad`: The largest gradient value for normalization.
- `para::parameters`: A struct containing various parameters, including the scale factor.

# Description
This function aligns the direction of an agent based on the chemotactic gradient. It performs the following steps:
1. Scales the agent's position based on the scale factor.
2. Ensures the scaled positions are not zero.
3. Calculates the gradient direction at the agent's position.
4. Computes the amplitude of the gradient and the angle between the agent's direction and the gradient direction.
5. Scales the angle by the amplitude of the gradient.
6. Constructs a rotation matrix based on the scaled angle.
7. Rotates the agent's direction using the rotation matrix.
8. Returns the new direction of the agent.

# Example
```julia
grad_x = rand(10, 10)
grad_y = rand(10, 10)
Y_pos = 5
X_pos = 5
X_dir = 1.0
Y_dir = 0.0
biggest_grad = 1.0
para = parameters(pa_ph = pa_ph(scale_fac = 1.0))

chemotactic_alignment(grad_x, grad_y, Y_pos, X_pos, X_dir, Y_dir, biggest_grad, para)
```
"""
function chemotactic_alignment(grad_x,grad_y,Y_pos,X_pos,X_dir,Y_dir,biggest_grad, para::parameters)
    y_x_scaled = round.(Int, [Y_pos, X_pos]./para.pa_ph.scale_fac .+0.01 )
    
    if y_x_scaled[1] == 0
        y_x_scaled[1] = 1
    end
    if y_x_scaled[2] == 0
        y_x_scaled[2] = 1
    end

    grad_dir = -[grad_y[y_x_scaled...], grad_x[y_x_scaled...]]
    amplitude_grad = sqrt(grad_dir[1]^2 + grad_dir[2]^2)/biggest_grad* para.pa_ph.Δt_walker_min  #para.pa_ph.Δt_walker_min is (mis)used as a scaling factor for the gradient (i know, horrible coding but otherwise my old data would not be importable)
    angle_grad = angle_between_vec([X_dir,Y_dir],grad_dir)
    
    amplitude_grad = (para.pa_ph.Δt_walker * amplitude_grad) > 1 ? 1 : para.pa_ph.Δt_walker * amplitude_grad
    scaled_angle = angle_grad*amplitude_grad
    rot_mat = [cos(scaled_angle) -sin(scaled_angle); sin(scaled_angle) cos(scaled_angle)]
    new_dir = rot_mat*[X_dir,Y_dir]
    return new_dir
end
