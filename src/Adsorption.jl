"""
    adsorb!(grids::Tuple{Matrix{Int64}, Matrix{Float64}}, full_agent_list, para::parameters, time_step::Int64 = 1)

Updates the gradient grid based on the positions of agents and their adsorption rates.

# Arguments
- `grids::Tuple{Matrix{Int64}, Matrix{Float64}}`: A tuple containing two matrices.
     The second matrix represents the concentration of a chemical.
- `full_agent_list`: A list of agents, which can be either a CircularBuffer of vectors of agents 
    or a vector of vectors of agents.
- `para::parameters`: A struct containing various parameters, including adsorption rate and 
    time step for the agents.
- `time_step::Int64`: The current time step (default is 1).

# Description
This function updates the gradient grid by increasing the concentration  at the positions of the agents.
     The amount of increase is determined by the adsorption rate and the time step specified in the `para` parameter.

"""
function adsorb!(grids::Tuple{Matrix{Int64}, Matrix{Float64}},full_agent_list,para::parameters, time_step::Int64 = 1)
    
    T = grids[2]

    if typeof(full_agent_list) == CircularBuffer{Vector{agent}}
        agent_list = full_agent_list[end]
    else
        agent_list = full_agent_list[time_step + 1]
    end

    for i = 1:length(agent_list)
        X_pos = agent_list[i].x_pos
        Y_pos = agent_list[i].y_pos

        y_x_scaled = round.(Int, [Y_pos, X_pos]./para.pa_ph.scale_fac .+0.01 )
    
        if y_x_scaled[1] == 0
            y_x_scaled[1] = 1
        end
        if y_x_scaled[2] == 0
            y_x_scaled[2] = 1
        end

        T[y_x_scaled...] += para.pa_ph.adsorption_rate * para.pa_ph.Δt_walker 
    end
end
