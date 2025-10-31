
"""
    create_para_struct(settings::Dict) -> parameters::parameters

Constructs a `parameters` struct from a dictionary of settings, dynamically generating the struct fields
     and their values based on the input dictionary. This function uses metaprogramming to parse and evaluate a 
     string representation of the struct initialization.

# Arguments
- `settings::Dict`: A dictionary where keys are the names of the parameters (as symbols or strings) and values are the
 corresponding values for these parameters. The values can be of any type, including `String`.

# Returns
- A new instance of the `parameters` struct with fields and values corresponding to the entries in the `settings` dictionary.
 If some keys are missing in the dictionary, the default values from the `parameters` struct definiton  will be used. 

# Examples
```julia
settings = Dict(:radius_collision => 2, :size => (10,))
para_struct = create_para_struct(settings)
```

This will create an instance of the parameters struct with radius_collision set to 2, size set to (10,10).

"""
function create_para_struct(settings::Dict)
    k = keys(settings)
    #create long defining string of Parameter struct
    para_ex = "parameters_physical("
    for i in k
        if typeof(settings[i]) != String
            para_ex = para_ex*"$i = $(settings[i]),"
        elseif typeof(settings[i]) == String
            para_ex = para_ex*"$i =\"$(settings[i])\","
        end
    end
    para_ex = chop(para_ex)*")"
    # parse the string with meta programming to create struct
    para_struct = eval(Meta.parse(para_ex))
    return para_struct
end

"""
    make_para_dict_list(sweep_dict_original::Dict)

Generate a list of parameter dictionaries for simulation runs based on the original sweep parameters dictionary.

# Arguments
- `sweep_dict_original::Dict`: A dictionary containing the original sweep parameters. 
It must include keys for `sweep_parameter`, `parameter_steps`,
 `iterations`, and `parameter_increment`. The `sweep_parameter` key should match one of the keys in the dictionary 
 that will be varied across simulations.

# Returns
- `Array{Dict}`: An array of dictionaries, each representing a set of parameters for a single simulation run.
 The specified sweep parameter is incremented across these dictionaries according to the `parameter_steps` 
 and `parameter_increment` specified in the original dictionary.

# Throws
- An error if the `sweep_parameter` specified in `sweep_dict_original` does not exist in the dictionary keys.

# Example
```julia
sweep_dict_original = Dict(
    "sweep_parameter" => "agent_number",
    "parameter_steps" => 5,
    "iterations" => 10,
    "parameter_increment" => 2,
    "other_param" => 1
)
para_dict_list = make_para_dict_list(sweep_dict_original)
```
This function is useful for creating parameter sweeps for simulations,
     where one wants to vary a single parameter across a range of values while keeping other parameters constant.  

"""
function make_para_dict_list(sweep_dict_original)
    sweep_dict = deepcopy(sweep_dict_original)
    if sweep_dict["sweep_parameter"] in keys(sweep_dict)
        para_dict_list = []
        para_delta = 0
        para_steps = pop!(sweep_dict,"parameter_steps")
        iterations = pop!(sweep_dict,"iterations")
        para_increment = pop!(sweep_dict,"parameter_increment")
        
        
        for i =1:para_steps
            for j = 1:iterations
                # copy dict per sim run without ensemble data
                para_temp = deepcopy(sweep_dict)
                delete!(para_temp,"sweep_parameter")
                # increase sweep para in run dict
                para_temp[sweep_dict["sweep_parameter"]] += para_delta
                push!(para_dict_list, para_temp)
            end
            para_delta += para_increment
        end            
    return para_dict_list
    end
   
    #else
       # throw("Used sweep_parameter key is not valid, use one of the following:boundary_conditions, agent_number, noise_dis,size, velocity,interaction_radius timesteps,noise_strength ")
   # end

end


"""
    remove_fields_dict(dict::Dict{String, Any}) -> Dict{String, Any}

Removes keys from the given dictionary `dict` that are not present in the `parameters` structure.

# Arguments
- `dict::Dict{String, Any}`: The dictionary from which keys will be removed.

# Returns
- `Dict{String, Any}`: The modified dictionary with only the keys that are present in the `parameters` structure.

# Example
```julia
parameters = (a=1, b=2)
dict = Dict("a" => 1, "b" => 2, "c" => 3)
new_dict = remove_fields_dict(dict)
# new_dict is now Dict("a" => 1, "b" => 2)
```
This function is useful for ensuring that a dictionary contains only valid keys before using it to set parameters in a simulation.
"""
function remove_fields_dict!(dict)
    fields = String.(fieldnames(parameters_physical))
    keys2 = keys(dict)  
    for key2 in keys2
        if key2 ∉ fields 
            pop!(dict, key2)
        end
    end
    return dict
end


"""
    sweep_and_save_parallel_t(para_dict_list, full_data_path, offset=100000; threads=4)

Run simulations in parallel and save the results.

# Arguments
- `para_dict_list::Vector{Dict}`: A list of parameter dictionaries for each simulation.
- `full_data_path::String`: The path where the simulation data will be saved.
- `offset::Int`: An offset for the file naming of the saved data. Default is 100000.
- `threads::Int`: The number of threads to use for parallel execution. Default is 4.

# Description
This function runs multiple simulations in parallel using the specified number of threads. 
    Each simulation is initialized with parameters from `para_dict_list`. 
    The simulation involves updating grids and agent lists over a series of time steps. 
    Data is sampled at specified intervals and stored in vectors. 
    After the simulation, the data is serialized and saved to `full_data_path` with filenames 
    based on the task index and offset.

# Example
```julia
para_dict_list = [...]  # List of parameter dictionaries
full_data_path = "path/to/save/data"
sweep_and_save_parallel_t(para_dict_list, full_data_path, offset=100000, threads=4) 
```
"""
function sweep_and_save_parallel_t(para_dict_list,full_data_path,offset =100000; threads = 4) 
    len = length(para_dict_list)
    Threads.@threads for o in 1:threads
        for j = o:threads:len
            println("initalize simulation")
            sampling_interval = pop!(para_dict_list[j],"sampling_interval")
            para = parameters(pa_ph = create_para_struct(remove_fields_dict!(deepcopy(para_dict_list[j]))))
            grids = create_grids(para)
            grids, Full_agent_list = initialize_system(grids,para)
            grids, Agent_list_cicular = make_sys_interactive!(grids,Full_agent_list)
            diff_grids = create_diff_grids(para)

            # create storage vec for agent and concentration grid
            vec_agent_storage = Vector{Vector{agent}}(undef,0)
            vec_gradient_grid_storage = Vector{Matrix{Float64}}(undef,0)
            vec_agent_grid_storage = Vector{Matrix{Int}}(undef,0)

            # add first entry at t=0
            push!(vec_agent_storage,deepcopy(Agent_list_cicular[end]))
            push!(vec_gradient_grid_storage,deepcopy(grids[2]))
            push!(vec_agent_grid_storage,deepcopy(grids[1]))   
            
            println("start simulation for ", para.pa_ph.total_time, " seconds which are: ", para.pa_ph.time_steps_to_compute, 
            " steps to compute on the gradient grid and \n",
            para.pa_ph.time_steps_to_compute/para.pa_ph.ratio_walker_diff ," agent steps to compute")

            t = 1
            sampling_counter = 0.0

            for i = 1:(para.pa_ph.time_steps_to_compute-1)
                diffusion_2D!(grids,diff_grids,para)
                
                if (i%para.pa_ph.ratio_walker_diff) == 0 
                    t += 1
                    

                    update_directions!(grids, Agent_list_cicular, para,t)
                    update_position!(grids, Agent_list_cicular, para,t)
                    strenghten_boundary!(grids, para)
                    divison!(grids,Agent_list_cicular,para,t)
                    adsorb!(grids, Agent_list_cicular, para,t)
                end
                if i%30 == 0
                    println("Outer Thread $o working on Task $j, stepcounter: ",i,"/",para.pa_ph.time_steps_to_compute)
                end

                sampling_counter += para.pa_ph.Δt_diff
                if sampling_counter >= sampling_interval
                    push!(vec_agent_storage,deepcopy(Agent_list_cicular[end]))
                    push!(vec_gradient_grid_storage,deepcopy(grids[2]))
                    push!(vec_agent_grid_storage,deepcopy(grids[1]))    
                    sampling_counter = 0.0
                end
            end
        

            
            data = (vec_agent_storage, vec_agent_grid_storage, vec_gradient_grid_storage, para)   
            
            
            s.serialize(full_data_path*"/"*string(j+offset)*".jls",data)

        end
    end
end



function sweep_and_save_parallel_t_static(para_dict_list,full_data_path,offset =100000; threads = 4) 
    len = length(para_dict_list)
    Threads.@threads for o in 1:threads
        for j = o:threads:len
            println("initalize simulation")
            sampling_interval = pop!(para_dict_list[j],"sampling_interval")
            para = parameters(pa_ph = create_para_struct(remove_fields_dict!(deepcopy(para_dict_list[j]))))
            grids = create_grids(para)
            grids, Agent_list_cicular, grid_vec_sparse = initialize_system_static(grids,para)
            diff_grids = create_diff_grids(para)

            # create storage vec for agent and concentration grid
            vec_agent_storage = Vector{Vector{agent}}(undef,0)
            vec_grid_storage = Vector{Matrix{Float64}}(undef,0)
            
            println("start simulation for ", para.pa_ph.total_time, " seconds which are: ", para.pa_ph.time_steps_to_compute, 
            " steps to compute on the gradient grid and \n",
            para.pa_ph.time_steps_to_compute/para.pa_ph.ratio_walker_diff ," agent steps to compute")

            t = 1
            sampling_counter = 0.0

            for i = 1:(para.pa_ph.time_steps_to_compute-1)
                diffusion_2D!(grids,diff_grids,para)
                
                if (i%para.pa_ph.ratio_walker_diff) == 0 
                    t += 1
                    grid_t = deepcopy(grids[1])

                    update_directions!(grids, Agent_list_cicular, para,t)
                    update_position!(grids, Agent_list_cicular, para,t)
                    update_grid!(grids,para, grid_t,grid_vec_sparse)
                    divison!(grids,Agent_list_cicular,para,t)
                    adsorb!(grids, Agent_list_cicular, para,t)
                end
                if i%30 == 0
                    println("Outer Thread $o working on Task $j, stepcounter: ",i,"/",para.pa_ph.time_steps_to_compute)
                end

                sampling_counter += para.pa_ph.Δt_diff
                if sampling_counter >= sampling_interval
                    push!(vec_agent_storage,deepcopy(Agent_list_cicular[end]))
                    push!(vec_grid_storage,deepcopy(grids[2]))
                    sampling_counter = 0.0
                end
            end
        

            
            data = (vec_agent_storage, grid_vec_sparse, vec_grid_storage, para)   
            
            
            s.serialize(full_data_path*"/"*string(j+offset)*".jls",data)

        end
    end
end