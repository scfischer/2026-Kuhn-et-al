function divison!(grids::Tuple{Matrix{Int64}, Matrix{Float64}},full_agent_list,para::parameters, timestep::Int = 1) 

    grid = grids[1]
    if typeof(full_agent_list) == CircularBuffer{Vector{agent}}
        agent_list = full_agent_list[end]
    else
        agent_list = full_agent_list[timestep+1]
    end
    number_agents = length(agent_list)

    for i in 1:length(agent_list)
        if rand() < para.growth_rate*para.pa_ph.Δt_walker
            sucess, X_pos_new, Y_pos_new = find_free_neighbour_points(grid, agent_list[i], para)
            if sucess
                Dir = rand(-π:π)
                number_agents += 1
                push!(agent_list, agent(X_pos_new, Y_pos_new, cos(Dir),sin(Dir)))

                grid[Y_pos_new,X_pos_new] = number_agents
            end
        end
    end
end



