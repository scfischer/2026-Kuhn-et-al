using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using Revise, Distributions, GLMakie
using TrypColonies
import Serialization as s

para_phys = parameters_physical(
    N                       = (1000,1000),
    L                       = (0.015,0.015),
    scale_fac               = 4,
    total_time              = 180.0,
    walker_speed_phy        = 5*10^-6,
    growth_rate             = 1.157*10^-5, 
    Diff_coe                = 7*10^-11,
    walker_step_size        = 3,
    adsorption_rate         = 0.05,
    decay_rate              = 0.01,
    Diameter_colony         = 0.003,
    agent_number            = 250000,
    noise_strength          = 1.2,
    grid_strength           = 300,  
    grid_recover_rate       = 30.0,
    radius_collision        = 12,
    radius_tanget           = 15,
    )

para = parameters( 
    pa_ph                   = para_phys,
    chemotaxis_flag         = true,
    repulsion_flag          = false,
    repulsion_range         = 2,
    geometry                = "circle",
    start_config_gradient   = "rand",
    )


    plot_para = plot_parameters(
    arrow_tip_size   = 6,
    arrow_tail_length= 0.00005,
    framerate        = 30, 
    fontsize         = 15,
    refresh_delay    = 0.005,
    res              = (900,1300),
    cor_func_scalar  = order_parameter,
    draw_arrows      = false,
    cor_func_vector  = radial_distribution,   
    type_cor_func    = "vector",
    timesteps_corfunc = 100,
    );

makie_config!(plot_para)
grids = create_grids(para)

grids, Full_agent_list, grid_vec_sparse = initialize_system(grids,para)
grids, Full_agent_list = make_sys_interactive!(grids, Full_agent_list)
diff_grids = create_diff_grids(para)

println("start simulation for ", para.pa_ph.total_time, " seconds which are: ", para.pa_ph.time_steps_to_compute, 
" steps to compute on the gradient grid and \n", para.pa_ph.time_steps_to_compute/para.pa_ph.ratio_walker_diff ," agent steps to compute")


#timestepping
let t = 1
    @time begin
        
        for i = 1:(para.pa_ph.time_steps_to_compute-1)
            diffusion_2D!(grids,diff_grids,para)
            
            if (i%para.pa_ph.ratio_walker_diff) == 0 
                t += 1
                grid_t = deepcopy(grids[1])

                update_directions!(grids, Full_agent_list, para,t)
                update_position!(grids, Full_agent_list, para,t)
                update_grid!(grids,para, grid_t,grid_vec_sparse)
                divison!(grids,Full_agent_list,para,t)
                adsorb!(grids, Full_agent_list, para,t)
            end
            if i%10 == 0
                println("stepcounter: ",i,"/",para.pa_ph.time_steps_to_compute)
            end
        end
    end
end


#save data

s.serialize("test.jls",(Full_agent_list,grid_vec_sparse,grids))
"""
f = Figure(size = (1400,1800))
display(f)
println("start Visualisation")

grids[1][:] = Matrix(grid_vec_sparse[1])

ax1,fig_obs = ini_animation(Full_agent_list,para,grids,f,plot_para);

@time begin
    record(f,"vids/static_test_diff5.mp4", 2:para.timesteps, framerate = 60) do i
        animation_step!(ax1,Full_agent_list,grid_vec_sparse, grids,fig_obs, plot_para, para)
    end
end
"""
    
