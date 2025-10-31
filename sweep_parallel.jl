using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using Revise, Distributions, GLMakie, CSV, Dates, DataFrames
using TrypColonies
import Serialization as s

doubling_time = 3600*9
growth_rate =  2^(1/doubling_time) -1

doubling_time_sec = 3600*20
growth_rate_sec = 2^(1/doubling_time_sec) -1
growth_rate_inc = growth_rate_sec - growth_rate

 settings = Dict{String,Any}(
    "iterations"            => 1,
    "parameter_steps"       => 2,
    "sweep_parameter"       => "growth_rate",
    "parameter_increment"   => growth_rate_inc,
    "sampling_interval"     => 300, # in seconds   

    "N"                     => (1000,1000),
    "scale_fac"             => 4,
    "L"                     => (0.015,0.015),
    "agent_number"          => 250000,
    "walker_speed_phy"      => 5*10^-6,
    "Diff_coe"              => 7*10^-11,
    "total_time"            => 37000.0,
    "growth_rate"           => growth_rate,
    "walker_step_size"      => 3,
    "decay_rate"            => 0.001,
    "adsorption_rate"       => 0.00001,
    "Diameter_colony"       => 0.003,

    "noise_dis"             => Normal,
    "noise_strength"        => 0.2,
    "velocity_dis"          => Normal,
    "velo_var_relative"     => 1.5,  

    "Δt_walker_min"         => 0.1, # turing ratio per second
    
    "grid_recover_rate"     => 7.0, ## corrected for \Delta t_walker 15/2.25
    "grid_strength"         => 8000,
    "radius_tanget"         => 15,
    "radius_collision"      => 16,

    )

#create timeseed for saving 
time_ms  = string(time())
time_sec = time_ms[end-8:end-5]


# choose if data should be stored locally in model folder or somewhere else on a static path
data_path = find_data_path(local_data = false)


name = "growth_rate_high_sample_rate"

#create folder, if it does not exist already
full_data_path = data_path*name*"_"*string(settings["agent_number"])*"_time_seed_"* time_sec 

mkpath(full_data_path)


# create settings and parameters txts, each one for manual reading and one jls file for import back into julia
para = parameters(pa_ph = create_para_struct(remove_fields_dict!(deepcopy(settings))))
s.serialize(full_data_path*"/parameters.jls",para)
open(full_data_path*"/parameters.txt","a") do io
    println(io,para)
 end
 
s.serialize(full_data_path*"/settings.jls",settings)
CSV.write(full_data_path*"\\settings.csv",settings)


# load database DataFrame


#df = DataFrame()

para_dict_list = make_para_dict_list(settings)

threads2 = round(Int,Threads.nthreads()/2)


t = @elapsed begin 
    sweep_and_save_parallel_t(para_dict_list,full_data_path,threads = threads2)
end
println("This run took $t seconds")

# modify and save settings in Database df
settings2 = deepcopy(settings)


settings2["path"] = full_data_path
settings2["date"] = unix2datetime(time())
settings2["unix_date"] = time() 
settings2["run_time"] = t
settings2["storage"] = foldersize(full_data_path)
settings2["location"] = data_path
settings2["name"] = name

df = s.deserialize("data_base/df.jls")   
println(last(df))

push!(df,settings2, promote = true, cols = :union )
println(last(df))
s.serialize("data_base/df.jls",df)
CSV.write("data_base/df.csv",df)

println(last(df))
println("finished")

