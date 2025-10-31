"""
    parameters_physical

A struct representing the physical parameters for a simulation involving walkers and diffusion processes. This struct is initialized with default values and includes various fields related to the simulation's spatial and temporal properties.

# Fields
- `N::Tuple{Int, Int}`: The number of grid points in the x and y directions in the gradient grid. Default is `(1000, 1000)`.
- `scale_fac::UInt`: The scaling factor between the gradient grid where the diffusion is calculated and the one of the ABM. Default is `4`.
- `L::Tuple{Float64, Float64}`: The physical dimensions of the simulation domain in the x and y directions. Default is `(0.015, 0.015)`.
- `agent_number::Int64`: The number of agents in the simulation. Default is `250000`.
- `walker_speed_phy::Float64`: The physical speed of the walkers in m/s. Default is `5*10^-6`.
- `Diff_coe::Float64`: The diffusion coefficient in m^2/s. Default is `7*10^-10`.
- `total_time::Float64`: The total simulation time in s. Default is `200.0`.
- `growth_rate::Float64`: The growth rate in 1/s. Default is `1.157*10^-5`.
- `walker_step_size::Int`: The step size for the walkers on the ABM grid. Default is `3`.
- `decay_rate::Float64`: The decay rate in 1/s. Default is `0.001`.
- `adsorption_rate::Float64`: The adsorption rate/amount in 1/s . Default is `0.001`.
- `Diameter_colony::Float64`: The diameter of the colony in m. Default is `0.003`.
- `noise_dis::Any`: The distribution for noise. Default is `Normal`.
- `noise_strength::Number`: The stdd of the noise distribution. Default is `0.5`.
- `velocity_dis::Any`: The distribution for velocity. Default is `Normal`.
- `velo_var_relative::Float64`: The relative variance in velocity. Default is `1.5`.
- `grid_strength::Int64`: The strength of the grid. Default is `100`.
- `grid_recover_rate::Float64`: The recovery rate of the grid. Default is `10.0`.
- `radius_tanget::Int64`: The radius to calculate the tangent which is needed to calculate the reflection on walls. Default is `6`.
- `radius_collision::Int64`: The radius for collision interactions. Default is `5`.
- `N_abm::Tuple{Int, Int}`: The number of grid points in the ABM grid, calculated as `N .* scale_fac`.
- `dx::Float64`: The grid spacing in the x direction, calculated as `L[1]/N[1]`.
- `dy::Float64`: The grid spacing in the y direction, calculated as `L[2]/N[2]`.
- `xc::LinRange{Float64, Int64}`: A linear range of x coordinates from `0.0` to `L[1]` with `N[1]` points.
- `yc::LinRange{Float64, Int64}`: A linear range of y coordinates from `0.0` to `L[2]` with `N[2]` points.
- `Δt_diff_min::Float64`: The minimum time step for diffusion, calculated as `(minimum((L./(N)))^2*2)/(8*Diff_coe+minimum((L./(N)))^2*decay_rate)` (see Von Neumann stability analysis).
- `Δt_walker_min::Float64`: The minimum time step for walkers, calculated as `(minimum((L./(N_abm)))*walker_step_size)/walker_speed_phy`.
- `Δt_diff::Float64`: The time step for diffusion, calculated using `scale_time_step`.
- `Δt_walker::Float64`: The time step for walkers, calculated using `scale_time_step`.
- `ratio_walker_diff::Int`: The ratio of walker time steps to diffusion time steps, calculated using `scale_time_step`.
- `time_steps_to_compute::Int`: The total number of time steps to compute, calculated as `round(Int, total_time/Δt_diff)`.
"""
@with_kw mutable struct parameters_physical
    N::Tuple{Int,Int}                   = (1000,1000)
    scale_fac::UInt                     = 4
    L::Tuple{Float64,Float64}           = (0.015,0.015)
    agent_number::Int64                 = 250000
    walker_speed_phy::Float64           = 5*10^-6
    Diff_coe::Float64                   = 7*10^-10
    total_time::Float64                 = 200.0
    
    growth_rate::Float64                = 1.9254273710078706e-5
    walker_step_size::Int               = 3
    decay_rate::Float64                 = 0.001
    adsorption_rate::Float64            = 0.001
    
    Diameter_colony::Float64            = 0.003
    
    noise_dis::Any                      = Normal 
    noise_strength::Number              = 0.5
    velocity_dis::Any                   = Normal
    velo_var_relative::Float64          = 1.5
    grid_strength::Int64                = 100
    grid_recover_rate::Float64          = 10.0
    radius_tanget::Int64                = 6
    radius_collision::Int64             = 5

    
    N_abm::Tuple{Int,Int}               = N .* scale_fac
    dx::Float64                         = L[1]/N[1]
    dy::Float64                         = L[2]/N[2]
    xc::LinRange{Float64, Int64}        = LinRange(0.0, L[1], N[1])
    yc::LinRange{Float64, Int64}        = LinRange(0.0, L[2], N[2])
    Δt_diff_min::Float64                = (minimum((L./(N)))^2*2)/(8*Diff_coe+minimum((L./(N)))^2*decay_rate)
    
    Δt_walker_min::Float64              = 1.0   
    Δt_walker::Float64                  = scale_time_step(Δt_diff_min, (minimum((L./(N_abm)))*walker_step_size)/walker_speed_phy)[2]
    Δt_diff::Float64                    = scale_time_step(Δt_diff_min, (minimum((L./(N_abm)))*walker_step_size)/walker_speed_phy)[1]
    ratio_walker_diff::Int              = scale_time_step(Δt_diff_min, (minimum((L./(N_abm)))*walker_step_size)/walker_speed_phy)[3]
    time_steps_to_compute::Int          = round(Int,total_time/Δt_diff)
end

"""
    parameters

Defines the parameters for the simulation with default values and possible options for certain fields.

# Fields
- `pa_ph::parameters_physical`: A struct containing physical parameters for the simulation.
- `interaction_radius::Int64`: The radius within which agents interact. Default is 7.
- `boundary_conditions::String`: The type of boundary conditions. Default is "periodic".
- `geometry::String`: The geometry of the simulation area. Default is "circle".
- `start_config::String`: The initial configuration of agents. Default is "random".
- `start_config_gradient::String`: The initial configuration of the gradient. Default is "random".
- `boundary_moveable::Bool`: Whether the boundary is moveable. Default is true.
- `path_tracing::Bool`: Whether to trace the path of agents. Default is true.
- `sliding_movement::Bool`: Whether the boundary allows sliding. Default is false.
- `repulsion_flag::Bool`: Whether repulsion is enabled. Default is false.
- `repulsion_range::Int64`: The range of repulsion. Default is 3.
- `chemotaxis_flag::Bool`: Whether chemotaxis is enabled. Default is true.
- `size::Tuple{Int, Int}`: The dimensions of the simulation grid. Default is `pa_ph.N_abm`.
- `xc_w::LinRange{Float64, Int64}`: A linear range of x coordinates for the grid.
- `yc_w::LinRange{Float64, Int64}`: A linear range of y coordinates for the grid.
- `arrow_to_grid_fac_x::Float64`: The scaling factor for x coordinates from arrows to grid.
- `arrow_to_grid_fac_y::Float64`: The scaling factor for y coordinates from arrows to grid.
- `velocity::Float64`: The initial velocity of agents. Default is `pa_ph.walker_step_size`.
- `timesteps::Int64`: The number of timesteps the simulation runs. Default is calculated from `pa_ph.total_time` and `pa_ph.Δt_walker`.
- `velocity_variance::Float64`: The variance in agent velocities. Default is `pa_ph.walker_step_size * pa_ph.velo_var_relative`.
- `growth_rate::Float64`: The growth rate. Default is `pa_ph.growth_rate`.

## Usage Notes
- Use `reflective` boundary conditions only when agents cannot leave the grid space; otherwise, the simulation may crash silently.
- The `geometry` field allows for simulating different environmental layouts, which can significantly affect agent behavior and interactions.
- The `start_config` field determines the initial placement of agents, which can influence the dynamics of the simulation from the outset.

This struct is used to configure and initialize the simulation environment, affecting how agents move, interact, and evolve over time.
"""
@with_kw mutable struct parameters
    pa_ph::parameters_physical      = parameters_physical()
    interaction_radius::Int64       = 7
    boundary_conditions::String     = "periodic"
    geometry::String                = "circle"
    start_config::String            = "random"
    start_config_gradient::String   = "random"

    boundary_moveable::Bool         = true 
    path_tracing::Bool              = true
    sliding_movement::Bool          = false
    repulsion_flag::Bool            = false 
    repulsion_range::Int64          = 3
    chemotaxis_flag::Bool           = true
    size::Tuple{Int,Int}            = pa_ph.N_abm
    xc_w::LinRange{Float64, Int64}  = LinRange(0.0, pa_ph.L[1], pa_ph.N[1]*pa_ph.scale_fac)
    yc_w::LinRange{Float64, Int64}  = LinRange(0.0, pa_ph.L[2], pa_ph.N[2]*pa_ph.scale_fac)
    arrow_to_grid_fac_x::Float64    = size[1]/pa_ph.L[1]
    arrow_to_grid_fac_y::Float64    = size[2]/pa_ph.L[2]
    velocity::Float64               = pa_ph.walker_step_size
    timesteps::Int64                = round(Int,pa_ph.total_time/(pa_ph.Δt_walker))
    velocity_variance::Float64      = pa_ph.walker_step_size * pa_ph.velo_var_relative
    growth_rate::Float64            = pa_ph.growth_rate
end



"""
    scale_time_step(dt_diff::Float64, dt_walker::Float64) -> Tuple{Float64, Float64}

Scales the time step `dt_diff` to be compatible with `dt_walker`.

# Arguments
- `dt_diff::Float64`: The differential time step.
- `dt_walker::Float64`: The walker time step.

# Returns
- `Tuple{Float64, Float64}`: A tuple containing the scaled `dt_diff` and the original `dt_walker`.

# Description
This function adjusts the differential time step `dt_diff` to ensure it is compatible with the walker time step `dt_walker`. If `dt_diff` is greater than or equal to `dt_walker`, it is set to `dt_walker`. Otherwise, it is adjusted to be a divisor of `dt_walker` with a tolerance for floating-point precision.
"""
function scale_time_step(dt_diff,dt_walker)
    if dt_diff >= dt_walker
        dt_diff = dt_walker
    else
        if dt_walker % dt_diff < 0.00000001 
            dt_diff = dt_walker/(dt_walker÷dt_diff)
        else
            dt_diff = dt_walker/(dt_walker÷dt_diff +1)
        end
    end
    return dt_diff, dt_walker, round(Int,dt_walker/dt_diff)
end


"""
    agent

A struct representing an agent in a grid. It holds the agent's current position in the grid as well as its direction of movement.

# Fields
- `x_pos::Int64`: The x-coordinate of the agent's position in the grid.
- `y_pos::Int64`: The y-coordinate of the agent's position in the grid.
- `dir_x::Float64`: The x-component of the agent's direction vector.
- `dir_y::Float64`: The y-component of the agent's direction vector.
"""
struct agent
    x_pos::Int64
    y_pos::Int64
    dir_x::Float64
    dir_y::Float64
end  

"""
    create_grid(para::parameters) -> Matrix{Int64}

Create a grid based on the specified parameters. The function initializes a grid of zeros and modifies it according to 
the geometry specified in `Parameters`. The grid can be shaped into a square, channel, circle, or a smaller circle,
with negative values (defined by `Parameters.grid_strength`) in the entire grid except for the specified geometry shape.

# Arguments
    - `para::parameters`: A struct containing the grid parameters. It must have the following fields:
        - `size::Tuple{Int64, Int64}`: The dimensions of the grid.
        - `geometry::String`: The shape of the grid, which can be "square", "channel", "circle", or "circle_small".
        - `grid_strength::Int64`: The negative value to be applied to the grid's border or outside the specified shape.

# Returns
- `Matrix{Int64}`: The generated grid with the specified geometry and dimensions.

# Examples
```julia
params = parameters(size=(100, 100), geometry="square", grid_strength=5)
grid = create_grid(params)
```
This function supports creating grids with different geometries to simulate various environments or constraints. 
    """ 
function create_grids(Parameters::parameters; create_scal_grid = true)
    size = Parameters.size
    grid = zeros(Int64,size[1],size[2])
    
    if Parameters.geometry =="square"
        border_size = size ÷ 10
        grid .= -Parameters.pa_ph.grid_strength
        grid[border_size+1:end-border_size,border_size+1:end-border_size] .= 0    
    
    elseif Parameters.geometry == "channel"
        border_size = size .÷ 10
        grid .= -Parameters.pa_ph.grid_strength
        grid[border_size[1]+1:end-border_size[1],1:end] .= 0   
    
    elseif Parameters.geometry == "circle"

        grid_strength = Parameters.pa_ph.grid_strength
        grid_deviation = round(Int, 0.2*grid_strength)
        
        grid -= rand((grid_strength-grid_deviation):(grid_strength+grid_deviation),size[1],size[2])
        #grid .= -Parameters.pa_ph.grid_strength
        create_circle!(grid,Parameters)
    end

    if create_scal_grid
        if Parameters.start_config_gradient == "exponential"
            scal_grid = create_exp_scalar_grid(Parameters)
        elseif Parameters.start_config_gradient == "binary"
            scal_grid = create_binary_scalar_grid(Parameters)
        else
            scal_grid = rand(Float64,Parameters.pa_ph.N[1],Parameters.pa_ph.N[2])*0.0002
            #scal_grid = zeros(Float64,Parameters.pa_ph.N[1],Parameters.pa_ph.N[2])
        end
    
        return (grid, scal_grid)
    else
        return grid
    end
end



function create_circle!(grid::Matrix{Int64},para::parameters)
    L = para.pa_ph.L
    diameter = para.pa_ph.Diameter_colony
    relative_diameter = diameter/L[1]

    length_grid = minimum(size(grid))
    center = [para.size[1] ÷ 2, para.size[2] ÷ 2]
    radius = round(Int64,(length_grid*relative_diameter)/2)

    for y = -radius:radius
        for x = -radius:radius
            if sqrt(x^2+y^2) <= radius
                grid[y+center[2],x+center[1]] = 0
            end
        end
    end
end



"""
    initialize_system!(grids::Tuple{Matrix{Int64}, Matrix{Float64}}, Parameters::parameters)

Initializes the simulation system by setting up the initial configuration of agents and grid states.

# Arguments
- `grids::Tuple{Matrix{Int64}, Matrix{Float64}}`: A tuple containing two matrices. The first matrix represents the grid to be initialized.
- `Parameters::parameters`: A struct containing various parameters, including the number of timesteps and other simulation parameters.

# Description
This function initializes the system by setting up the initial configuration of agents and grid states. It performs the following steps:
1. Initializes an empty list of agents.
2. Creates a list to store the state of agents for each timestep.
3. Calls `random_start_config!` to set up the initial configuration of the grid and agents.
4. Copies the initial agent list to all timesteps.
5. Initializes a sparse matrix vector to save grid changes.
6. Saves the initial grid state to the sparse matrix vector.

# Example
```julia
grids = (Matrix{Int64}(undef, 10, 10), Matrix{Float64}(undef, 10, 10))
Parameters = parameters(pa_ph = pa_ph(time_steps_to_compute = 100), timesteps = 100)

initialize_system(grids, Parameters)
```
# Returns
grids: The initialized grids.
agent_list_full: A list of agent states for each timestep.
grid_vec_sparse: A vector of sparse matrices to save grid changes. 
"""	
function initialize_system(grids::Tuple{Matrix{Int64}, Matrix{Float64}},Parameters::parameters)
    ini_length = 10
    # initalize configuration in first timestep
    agent_list  = Vector{agent}(undef, 0)
    #agent_list_full = Vector{Vector{agent}}(undef, Parameters.pa_ph.time_steps_to_compute)
    agent_list_full = Vector{Vector{agent}}(undef, ini_length)
    
    random_start_config!(grids[1],Parameters,agent_list)
    
    # initialize full agent list with copies of system state in first timestep
    for i = 1:ini_length
        agent_list_full[i] = deepcopy(agent_list)
    end

    ##  grid vector parse to save grid changes, only gets used in non interactive simulations
    grid_vec_sparse = Vector{SparseMatrixCSC{Int64, Int64}}(undef,0)

    #Initialize first entry of sparse grid vector with grid state at t=0    
    save_grid_changes!(zeros(Int, size(grids[1])),replace( u -> (u>0) ? 0 : u, grids[1]),grid_vec_sparse)
    
    return grids, agent_list_full, grid_vec_sparse
end

function initialize_system_static(grids::Tuple{Matrix{Int64}, Matrix{Float64}},Parameters::parameters)
    #ini_length = 10
    # initalize configuration in first timestep
    agent_list  = Vector{agent}(undef, 0)
    agent_list_full = Vector{Vector{agent}}(undef, Parameters.pa_ph.time_steps_to_compute)
    #agent_list_full = Vector{Vector{agent}}(undef, ini_length)
    
    random_start_config!(grids[1],Parameters,agent_list)
    
    # initialize full agent list with copies of system state in first timestep
    for i = 1:Parameters.pa_ph.time_steps_to_compute
        agent_list_full[i] = deepcopy(agent_list)
    end

    ##  grid vector parse to save grid changes, only gets used in non interactive simulations
    grid_vec_sparse = Vector{SparseMatrixCSC{Int64, Int64}}(undef,0)

    #Initialize first entry of sparse grid vector with grid state at t=0    
    save_grid_changes!(zeros(Int, size(grids[1])),replace( u -> (u>0) ? 0 : u, grids[1]),grid_vec_sparse)
    
    return grids, agent_list_full, grid_vec_sparse
end

"""
### `make_sys_interactive!`

Transforms a grid system into an interactive one by initializing a circular buffer to store agent states, allowing for dynamic updates.

#### Parameters:
- `grids`: The grids on which agents operate. This function does not modify the grid but returns it for consistency in the interface.
- `agent_list_full`: A list (vector) of agents. These agents are the initial state of the system.
- `buffer_length::Int=5` (optional): The length of the circular buffer. This determines how many past states of the agent list can be stored. Default is 5.

#### Returns:
- `grids`: The same grid passed as input, returned unchanged.
- `agent_list_cicular`: A circular buffer filled with deep copies of the first agent in `agent_list_full`.
 This buffer allows for storing and accessing past states of the system.

"""
function make_sys_interactive!(grids::Tuple{Matrix{Int64}, Matrix{Float64}}, agent_list_full; buffer_length = 5)
    #creat circular buffer
    agent_list_cicular = CircularBuffer{Vector{agent}}(buffer_length)

    for i in 1:buffer_length
        push!(agent_list_cicular,deepcopy(agent_list_full[1]))
    end

    return grids, agent_list_cicular
end

function random_start_config!(grid::AbstractArray{Int}, Parameters::parameters, agent_list::Vector{agent})
    agent_number = Parameters.pa_ph.agent_number
    x_size = size(grid)[2]
    y_size = size(grid)[1]
    i = 1
    while i <= agent_number
        X_pos = rand(1:x_size)
        Y_pos = rand(1:y_size)
        Dir = rand(-π:.001:π)
        Dir_x = cos(Dir)
        Dir_y = sin(Dir)
        if grid[Y_pos,X_pos] != 0
            continue
        end
        grid[Y_pos,X_pos] = i
        push!(agent_list,agent(X_pos,Y_pos,Dir_x,Dir_y))
        i = i+1
    end
end


"""
    reconstruct_grid!(Grid_vec_sparse::Vector{SparseMatrixCSC{Int64, Int64}}, Timestep::Int, start_grid::AbstractArray{Int})

Update `start_grid` in-place by adding the non-zero values from a specific timestep's sparse matrix in `Grid_vec_sparse`.

# Arguments
- `Grid_vec_sparse`: A vector of sparse matrices (`SparseMatrixCSC{Int64, Int64}`), each representing the grid at a different timestep.
- `Timestep`: The specific timestep (index into `Grid_vec_sparse`) from which to take the non-zero values.
- `start_grid`: An array to be updated in-place. It represents the grid before updating, and it will be modified to include the non-zero values from the specified timestep's sparse matrix.

# Description
This function iterates over all non-zero elements of the sparse matrix at the specified `Timestep` in `Grid_vec_sparse`. For each non-zero element, it increments the corresponding element in `start_grid` by the value of the non-zero element. This operation modifies `start_grid` in-place, effectively reconstructing or updating the grid state by applying the changes from the specified timestep.

# Examples
```julia
# Example usage
Grid_vec_sparse = [sparse([1, 2], [1, 2], [1, 1]), sparse([1, 2], [1, 2], [2, 2])]
Timestep = 2
start_grid = zeros(Int, 2, 2)

reconstruct_grid!(Grid_vec_sparse, Timestep, start_grid)
# start_grid is now [0 2; 0 2] (2x2 matrix of integers
```
"""
function reconstruct_grid!(Grid_vec_sparse::Vector{SparseMatrixCSC{Int64, Int64}},Timestep::Int, start_grid::AbstractArray{Int})
    non_zeros_index = findnz(Grid_vec_sparse[Timestep])
    for i in 1:length(non_zeros_index[1])
        start_grid[non_zeros_index[1][i],non_zeros_index[2][i]] += non_zeros_index[3][i]
    end
end


"""
    reconstruct_grid_from_scratch(Grid_vec_sparse::Vector{SparseMatrixCSC{Int64, Int64}}, Timestep::Int = length(Grid_vec_sparse)) -> Matrix{Int64}

Reconstructs the grid state from scratch up to a specified timestep using a sequence of sparse matrices.

# Arguments
- `Grid_vec_sparse`: A vector of `SparseMatrixCSC{Int64, Int64}` objects, each representing the grid changes at a specific timestep.
- `Timestep`: (Optional) The timestep up to which the grid should be reconstructed. Defaults to the last available timestep in `Grid_vec_sparse`.

# Returns
A `Matrix{Int64}` representing the reconstructed grid state up to the specified timestep.

# Description
This function reconstructs the grid state by starting with the grid state at the first timestep
     and then applying changes from each subsequent timestep up to the specified `Timestep`. 
     It uses the `reconstruct_grid!` function to apply changes in-place to the grid state.

"""
function reconstruct_grid_from_scratch(Grid_vec_sparse::Vector{SparseMatrixCSC{Int64, Int64}},Timestep::Int = length(Grid_vec_sparse))
    start_grid = Matrix(Grid_vec_sparse[1])
    for i = 2:Timestep
        reconstruct_grid!(Grid_vec_sparse,i,start_grid)
    end
    return start_grid
end


"""
    save_grid_changes!(Grid_1::Union{Matrix{Int64}, Real}, Grid_2::Matrix{Int64}, Grid_vec::Vector{SparseMatrixCSC{Int64, Int64}})

Subtracts `Grid_1` from `Grid_2` and appends the result to `Grid_vec`.

# Arguments
- `Grid_1`: A matrix of integers or a real number. Represents the initial state of the grid.
- `Grid_2`: A matrix of integers. Represents the updated state of the grid.
- `Grid_vec`: A vector of sparse integer matrices. This vector is updated with the changes.

# Notes
- This function modifies `Grid_vec` in place by appending the difference between `Grid_2` and `Grid_1`.
- If `Grid_1` is a real number, it is subtracted from each element of `Grid_2` before the result is appended to `Grid_vec`.
"""
function save_grid_changes!(Grid_1::Matrix{Int64},Grid_2::Matrix{Int64}, Grid_vec::Vector{SparseMatrixCSC{Int64, Int64}})
    grid_1 = replace(u -> (u>0) ? 0 : u, Grid_1)
    grid_2 = replace(u -> (u>0) ? 0 : u, Grid_2)
    push!(Grid_vec,grid_2 .-grid_1)
end

"""
    strenghten_boundary!(grids::Tuple{Matrix{Int64}, Matrix{Float64}}, Parameters::parameters)

Strengthens the boundary of the grid by recovering the grid values based on a recovery rate.

# Arguments
- `grids::Tuple{Matrix{Int64}, Matrix{Float64}}`: A tuple containing two matrices. The first matrix represents the grid to be strengthened.
- `Parameters::parameters`: A struct containing various parameters, including grid size, recovery rate, and grid strength.

# Description
This function iterates over the grid and increases the grid values that are negative but greater than the negative grid strength.
     The increase is determined by a Poisson distribution with a mean equal to the grid recovery rate. 
     If the updated grid value exceeds the negative grid strength, it is capped at the negative grid strength.

# Example
```julia
grids = (Matrix{Int64}(undef, 10, 10), Matrix{Float64}(undef, 10, 10))
Parameters = parameters(size = (10, 10), grid_recover_rate = 0.1, grid_strength = 5)

strenghten_boundary!(grids, Parameters)
```
"""
function strenghten_boundary!(grids::Tuple{Matrix{Int64}, Matrix{Float64}}, Parameters::parameters)
    grid = grids[1] 
    size = Parameters.size
    if Parameters.pa_ph.grid_recover_rate > 0.0 
        
        for y in 1:size[1]
            for x in 1:size[2]
                if grid[y,x] < 0 && grid[y,x] > -Parameters.pa_ph.grid_strength
                    grid_add = rand(Poisson(Parameters.pa_ph.grid_recover_rate * Parameters.pa_ph.Δt_walker))
                    if (grid[y,x] -= grid_add) < -Parameters.pa_ph.grid_strength
                        grid[y,x] = -Parameters.pa_ph.grid_strength
                    else
                        grid[y,x] -= grid_add
                    end
                end
            end
        end
    end
end


"""
    update_grid!(grid::Matrix{Int}, Parameters::parameters, grid_t::Matrix{Int}, grid_vec_sparse::Vector{SparseMatrixCSC{Int64, Int64}})

Updates the simulation grid by strengthening boundaries and saving grid changes.

# Arguments
- `grid`: The current state of the simulation grid as a matrix of integers.
- `Parameters`: A struct containing simulation parameters, such as grid size and boundary strength.
- `grid_t`: A reference grid state used for comparison when saving changes.
- `grid_vec_sparse`: A vector of sparse matrices where changes to the grid are stored.

# Description
This function performs two main operations on the simulation grid:
1. It calls `strengthen_boundary!` to modify the grid's boundaries based on the `Parameters`. This typically involves making the boundaries more negative (stronger) according to a predefined strength parameter.
2. It calls `save_grid_changes!` to compare the current grid state (`grid`) with a reference state (`grid_t`) and records the differences in `grid_vec_sparse`. This step is crucial for tracking the evolution of the grid over time in an efficient manner.
"""
function update_grid!(grids::Tuple{Matrix{Int64}, Matrix{Float64}}, Parameters::parameters, grid_t::Matrix{Int}, grid_vec_sparse::Vector{SparseMatrixCSC{Int64, Int64}})
    #update grid
    strenghten_boundary!(grids, Parameters)
    grid = grids[1]
    save_grid_changes!(grid_t,grid,grid_vec_sparse)
end



"""
    create_exp_scalar_grid(para::parameters) -> Matrix{Float64}

Creates a scalar grid based on a Gaussian distribution centered in the middle of the grid. The gradient values are calculated using the specified parameters.

# Arguments
- `para::parameters`: A struct containing the parameters for the grid, including:
  - `size::Tuple{Int, Int}`: The dimensions of the grid.

# Returns
- `Matrix{Float64}`: A 2D array representing the gradient grid, with values following a Gaussian distribution.

# Example Usage
```julia
parameters = parameters(size=(100, 100))
grad_grid = create_exp_scalar_grid(parameters)
```
This function is useful for creating a scalar field that can be used in
simulations where agents respond to spatial gradients, such as chemotaxis or other gradient-based behaviors. 
""" 
function create_exp_scalar_grid(para::parameters; central_factor = 0.5)
    N_x = para.pa_ph.N[1]
    N_y = para.pa_ph.N[2]
    con_fac =  N_x/5 

    xc = LinRange(0,N_x,N_x)
    yc = LinRange(0,N_y,N_y)

    scal_grid = exp.(.-((xc .- N_y*central_factor)/con_fac).^2 .-((yc' .- N_x/2)/con_fac).^2)
    return scal_grid
end



function create_binary_scalar_grid(para::parameters)
    grid = zeros(Float64, para.pa_ph.N[1], para.pa_ph.N[2])
    L = para.pa_ph.L
    diameter = para.pa_ph.Diameter_colony
    relative_diameter = diameter/L[1]

    length_grid = minimum(size(grid))
    center = [para.size[1] ÷ 2, para.size[2] ÷ 2]
    radius = round(Int64,(length_grid*relative_diameter)/2)

    for y = -radius:radius
        for x = -radius:radius
            if sqrt(x^2+y^2) <= radius
                grid[y+center[2],x+center[1]] = 2.0
            end
        end
    end
    return grid
end
