"""
    plot_parameters

A structure defining parameters for plotting with default values.

# Fields
- `framerate`: The number of frames per second in the animation. Default is 60.
- `theme`: The plotting theme, using `Attributes`. Default is `theme_ggplot2()`.
- `fontsize`: The font size for text elements. Default is 22.
- `arrow_tip_size`: The size of the arrow tips. Default is 12.
- `arrow_tail_length`: The length of the arrow tails. Default is 5.
- `refresh_delay`: The delay between refreshes in seconds. Default is 0.05.
- `res`: The resolution of the plot as a tuple (width, height). Default is (1400, 2000).
- `Color_arrows`: The color scheme for arrows. Default is the full range of `ColorSchemes.hsv`.
- `Color_heatmap`: The color scheme for heatmaps. Default is indices 50 to 100 of `ColorSchemes.bone`.
- `timesteps_corfunc`: The number of timesteps for the correlation function. Default is 100.
- `cor_func`: The correlation function for circles. Default is `order_parameter_circle`.

# Usage
This struct is used to configure the visual aspects and functionalities of plots, including animations, color schemes, and correlation calculations.
"""
@with_kw mutable struct plot_parameters
    framerate::Int64                            = 60
    theme::Attributes                           = theme_ggplot2()
    fontsize::Int64                             = 22
    arrow_tip_size::Float64                     = 12
    arrow_tail_length::Float64                  = 5 
    refresh_delay::Float64                      = 0.005
    res::Tuple{Int64,Int64}                     = (1400,2000)
    Color_arrows::Any                           = ColorSchemes.hsv[1:end]
    Color_heatmap::Any                          = ColorSchemes.bone[50:100]
    timesteps_corfunc::Int64                    = 100
    type_cor_func::String                       = "scalar" 
    cor_func_scalar::Function                   = order_parameter
    cor_func_vector::Function                   = radial_distribution
    alpha_heatmap::Float64                      = 0.5
    draw_arrows::Bool                           = true
    possible_cor_func::Vector{Vector{Any}}      = [[order_parameter, "scalar"], [order_parameter_circle, "scalar"],
                                                [radial_distribution, "vector"], [radial_density_c, "vector"],
                                                [angular_metric, "vector"],[area_gain, "scalar"],
                                                [adsorbed_material, "scalar"],[order_parameter_perpendic, "scalar"],[pair_cor_metric, "vector"]]
                                                
                                            
end

## Dictionary used for the plot legend
fu_dict = Dict(order_parameter => "v_a[%]", order_parameter_circle => "v_c[%]",
        radial_distribution => "radial density",radial_density_c =>"radial_density_c", angular_metric => "angular metric", area_gain => "area gain",
         adsorbed_material => "adsorbed mat",
        order_parameter_perpendic => "v_perp[%]", pair_cor_metric => "pair correlation")  

#type_cor_func = "scalar" or "vector" 
 

function makie_config!(plot_parameters::plot_parameters )
    set_theme!(plot_parameters.theme)
    fontsize_theme = Theme(fontsize = plot_parameters.fontsize)
    update_theme!(fontsize_theme)
    GLMakie.activate!()
    Makie.inline!(false)
end

"""This function counteracts the dynamic scaling of colormaps in Makie if min and max values change over the course of an animation.
 It is given a Max value (min value = 0) and shrinks the colormap dynamically if the data does not reach the max value. """
function scale_color_map(colormap,data_array;Zero_offset=1*π, Max = 2*π)
    len_color_map = length(colormap)
    data_array = (data_array .+Zero_offset)./Max 
    
    min_percentage = minimum(data_array)
    max_percentage = maximum(data_array)
    
    min_index = round(Int64,min_percentage*len_color_map)
    if min_index == 0
        min_index = 1
    end
    max_index = round(Int64,max_percentage*len_color_map)
    
    return colormap[min_index:max_index]
end

"""
    mutable struct Figure_observables

A structure to hold various observables for visualizing simulation data.

# Fields
- `Time::Observable{Int64}`: An observable for the current time step.
- `V_a::Observable{CircularBuffer{Point2f}}`: An observable for a circular buffer of 2D points represent a correlation function.
- `V_a_length::Observable{Int64}`: An observable for the length of `V_a`.
- `X::Observable{Vector{Float64}}`: An observable for the x-coordinates of agents.
- `Y::Observable{Vector{Float64}}`: An observable for the y-coordinates of agents.
- `U::Observable{Vector{Float64}}`: An observable for the x-components of agent velocities.
- `V::Observable{Vector{Float64}}`: An observable for the y-components of agent velocities.
- `N::Observable{Int64}`: An observable for the number of agents.
- `Grid_contour::Observable{Matrix{Int64}}`: An observable for the grid which the agents move on.
- `Theta::Observable{Vector{Float64}}`: An observable for the angles of agents.
- `Isrunning::Observable{Bool}`: An observable to indicate if the simulation is running.
- `Color_scaled`: An observable for the scaled color values to make non-dynamic or dynamic colorscaling possible.
"""
mutable struct Figure_observables
    Time::Observable{Int64}
    Time_physical::Observable{Float64}
    V_a::Observable{CircularBuffer{Point2f}}
    V_a_length::Observable{Int64}
    X::Observable{Vector{Float64}}
    Y::Observable{Vector{Float64}}
    U::Observable{Vector{Float64}}
    V::Observable{Vector{Float64}}
    N::Observable{Int64}
    Grid_contour::Observable{Matrix{Int64}}
    Grid_gradient::Observable{Matrix{Float64}}
    Theta::Observable{Vector{Float64}}
    Isrunning::Observable{Bool}
    Color_scaled::Observable{Any}
end

function create_Figure_observables(Full_Agent_list,Grids::Tuple{Matrix{Int64}, Matrix{Float64}}
    ,para::parameters,Plot_parameter::plot_parameters)
    if Plot_parameter.type_cor_func == "scalar"
        v_a_list = CircularBuffer{Point2f}(Plot_parameter.timesteps_corfunc)
        push!(v_a_list, Point2f(1.0,(Plot_parameter.cor_func_scalar(Full_Agent_list[1],Grids,para))))
    elseif Plot_parameter.type_cor_func == "vector"
        dist_vec, dist_vec2 =  Plot_parameter.cor_func_vector(Full_Agent_list[1], Grids,para)
        v_a_list = CircularBuffer{Point2f}(length(dist_vec))
        append!(v_a_list, Point2f.(1:length(dist_vec),dist_vec))
    end
    
    time = Observable(1)
    time_physical = Observable(0.0)
    v_a = Observable(v_a_list) # make it an observable
    x = Observable([i.x_pos/para.arrow_to_grid_fac_x for i in Full_Agent_list[1]])
    y = Observable([i.y_pos/para.arrow_to_grid_fac_y for i in Full_Agent_list[1]])
    u = Observable([i.dir_x for i in Full_Agent_list[1]])
    v = Observable([i.dir_y for i in Full_Agent_list[1]])
    N = @lift(length($v))
    grid_contour = Observable(replace(u -> (u>0) ? 0 : u, Grids[1]))
    grid_gradient = Observable(Grids[2])
    isrunning = Observable(true) 
    theta = @lift((atan.($v ,$u )))
    color_scaled = @lift scale_color_map(Plot_parameter.Color_arrows,$theta)
    Fig_obs = Figure_observables(time,
                                time_physical,
                                v_a,
                                Plot_parameter.timesteps_corfunc,
                                x,
                                y,
                                u,
                                v, 
                                N,
                                grid_contour,
                                grid_gradient,
                                theta,
                                isrunning,
                                color_scaled
                                )
    return Fig_obs
end


function return_core_func(Plot_parameter::plot_parameters)
    if Plot_parameter.type_cor_func == "scalar"
        return Plot_parameter.cor_func_scalar
    elseif Plot_parameter.type_cor_func == "vector"
        return Plot_parameter.cor_func_vector
    end
end

"""
    create_buttons(f::Figure, fig_obs::Figure_observables, Parameters::parameters) -> (Button, Observable{parameters})

Creates interactive UI elements, including a start/stop button and a grid of sliders, to control simulation parameters in a graphical figure.

# Arguments
- `f::Figure`: The main figure where the UI elements will be added.
- `fig_obs::Figure_observables`: An observable object that includes the running state of the simulation.
- `Parameters::parameters`: A struct containing the initial parameters for the simulation.

# Returns
- `(Button, Observable{parameters})`: A tuple containing the start/stop button and an observable parameter struct that updates based on the slider values.

# Example Usage
```julia
f = Figure()
fig_obs = Figure_observables(Isrunning = Observable(false))
Parameters = parameters(...)  # Assume this is the parameters struct with initial values
pause_button, para_obs = create_buttons(f, fig_obs, Parameters)
```
"""
function create_buttons(f::Figure,fig_obs::Figure_observables,Parameters::parameters, Plot_parameter::plot_parameters)

    # Proxy grid for enought space for the buttons  
    f[3,:] = GridLayout(tellwidth = false, height = Plot_parameter.res[2]*0.02)

    f[4,:] = buttongrid = GridLayout(tellwidth = false)

   pause_button = Button(buttongrid[1,2]; label = " ⏯︎", tellwidth = false )

   
   menu = Menu(buttongrid[2,2], options =  zip([x[1] for x in Plot_parameter.possible_cor_func],Plot_parameter.possible_cor_func)
   , default = string(return_core_func(Plot_parameter)),  width = Plot_parameter.res[1]*0.2)

   on(menu.selection) do s
        type_cor_func = s[2]
        if type_cor_func == "scalar"
            Plot_parameter.cor_func_scalar = s[1]
            Plot_parameter.type_cor_func =  type_cor_func
            f.content[1].title = fu_dict[Plot_parameter.cor_func_scalar]
            #fig_obs.V_a[] = empty!(fig_obs.V_a[])
            fig_obs.V_a[] = CircularBuffer{Point2f}(Plot_parameter.timesteps_corfunc)

        elseif type_cor_func == "vector"
            Plot_parameter.cor_func_vector = s[1]
            Plot_parameter.type_cor_func = type_cor_func
            f.content[1].title = fu_dict[Plot_parameter.cor_func_vector]
            fig_obs.V_a[] = empty!(fig_obs.V_a[])
        end
    end

   pause_listener = on(pause_button.clicks) do shit 
       fig_obs.Isrunning[]= !fig_obs.Isrunning[]
       println("button clicked")
       println(fig_obs.Isrunning[])
       notify(fig_obs.Isrunning)
   end 


    lsgrid = SliderGrid(buttongrid[:,1],
    (label ="Noise [rad]", range =0.0:0.05:2.5,startvalue = Parameters.pa_ph.noise_strength),
    (label ="Interaction_range [N]", range =1:50,startvalue = Parameters.interaction_radius),
    (label ="Velocity [m/s]", range =0.0000001:0.0000001:0.00005,startvalue = Parameters.pa_ph.walker_speed_phy),
    (label ="Velocity_variance [%]", range =0:0.01:5.0,startvalue = Parameters.pa_ph.velo_var_relative),
    (label ="Diffusion coefficient [m²/s]",range  = collect(10 .^(range(-11,-8, length=400))), startvalue = Parameters.pa_ph.Diff_coe),
    (label ="Adsorption rate [1/s]", range =pushfirst!(collect(10 .^(range(-4,-1, length=400))), 0.0),startvalue = Parameters.pa_ph.adsorption_rate),
    (label = "growth_rate [1/s]", range = 0.0:0.000001:0.0001, startvalue = Parameters.growth_rate),
    (label = "decay_rate [1/s]", range = pushfirst!(collect(10 .^(range(-4,-1, length=400))), 0.0), startvalue = Parameters.pa_ph.decay_rate )
    , width = Plot_parameter.res[1]*0.6, tellwidth = false, halign = :left)

    #create new parameter struct including the observables that can change
    noise_obs = lsgrid.sliders[1].value
    interaction_radius_obs = lsgrid.sliders[2].value
    velocity_obs = lsgrid.sliders[3].value
    velocity_variance_obs = lsgrid.sliders[4].value
    diffusion_obs = lsgrid.sliders[5].value
    adsorption_rate_obs = lsgrid.sliders[6].value
    growth_rate_obs = lsgrid.sliders[7].value
    decay_rate_obs = lsgrid.sliders[8].value


   #create and link start/stop button to isrunning observable
   
   para_obs_physical = @lift(parameters_physical(
                                N                                   = Parameters.pa_ph.N,
                                walker_speed_phy                    = $velocity_obs,
                                L                                   = Parameters.pa_ph.L,
                                scale_fac                           = Parameters.pa_ph.scale_fac,
                                Diff_coe    	                    = $diffusion_obs,
                                total_time                          = Parameters.pa_ph.total_time,
                                growth_rate                         = $growth_rate_obs, 
                                walker_step_size                    = Parameters.pa_ph.walker_step_size,
                                decay_rate                          = $decay_rate_obs,
                                adsorption_rate                     = $adsorption_rate_obs,
                                Diameter_colony                     = Parameters.pa_ph.Diameter_colony,
                                agent_number                        = Parameters.pa_ph.agent_number,
                                velocity_dis                        = Parameters.pa_ph.velocity_dis,
                                velo_var_relative                   = $velocity_variance_obs,
                                grid_strength                       = Parameters.pa_ph.grid_strength,
                                grid_recover_rate                   = Parameters.pa_ph.grid_recover_rate,
                                radius_tanget                       = Parameters.pa_ph.radius_tanget,
                                radius_collision                    = Parameters.pa_ph.radius_collision,
                                noise_strength                      = $noise_obs,
                                ))

    para_obs = @lift(parameters(size                                = Parameters.size,
                                pa_ph                               = $para_obs_physical,
                                
                                interaction_radius                  = $interaction_radius_obs,
                                
                                timesteps                           = Parameters.timesteps,
                                boundary_conditions                 = Parameters.boundary_conditions,
                                geometry                            = Parameters.geometry,
                                
                                start_config                        = Parameters.start_config,
                                
                                
                                boundary_moveable                   = Parameters.boundary_moveable,
                                path_tracing                        = Parameters.path_tracing,
                                sliding_movement                    = Parameters.sliding_movement,
                                
                                
                                repulsion_flag                      = Parameters.repulsion_flag,
                                repulsion_range                     = Parameters.repulsion_range,
                                chemotaxis_flag                     = Parameters.chemotaxis_flag,
                                growth_rate                         = $growth_rate_obs,
                                ))

    return pause_button,para_obs

end



function ini_animation(Full_Agent_list,Parameters::parameters,Grids::Tuple{Matrix{Int64}, Matrix{Float64}}
    , f::Figure, Plot_parameter::plot_parameters)

    #create stuct with all Observable
    fig_obs = create_Figure_observables(Full_Agent_list,Grids,Parameters,Plot_parameter)

    #top lineplot axis 
    ax1 = Axis(f[1,1], title = "$(Plot_parameter.type_cor_func == "scalar" ? fu_dict[Plot_parameter.cor_func_scalar] : fu_dict[Plot_parameter.cor_func_vector])")
    rowsize!(f.layout, 1, Relative(0.15))
    line_1 = lines!(ax1,fig_obs.V_a)

    
    xlims!(ax1, fig_obs.Time[],fig_obs.V_a_length[]) 
    ylims!(ax1, -0.1, 1.1) 
    # main plot vector plot of agents: 
    ax2 = Axis(f[2, 1],title = @lift("Time[s] = $(round(($(fig_obs.Time_physical)), digits = 3)), N = $(fig_obs.N[]) "),
        xlabel = "X[m]",ylabel = "Y[m]",xlabelsize=20,ylabelsize=20)
    rowsize!(f.layout, 2, Relative(0.60))
    xlims!(ax2, 0, Parameters.pa_ph.L[2]) 
    ylims!(ax2, 0, Parameters.pa_ph.L[1])
    
    # create contour map out of grid which includes some kind of walls 
    if Parameters.geometry != "no_walls"

        ##change rotation to inversion and rotation
        heatmap!(ax2,Parameters.xc_w,Parameters.yc_w,fig_obs.Grid_contour,colormap = RGBAf.(Plot_parameter.Color_heatmap,Makie.alpha.(Plot_parameter.Color_heatmap) .* Plot_parameter.alpha_heatmap))
    end
    # add vector plot
    if Plot_parameter.draw_arrows == true
        arrows2d!(ax2,fig_obs.Y, fig_obs.X, fig_obs.V, fig_obs.U, 
        #arrowsize = Plot_parameter.arrow_tip_size, 

        lengthscale = Plot_parameter.arrow_tail_length,
        #arrowcolor = fig_obs.Theta, linecolor = fig_obs.Theta,
        color = fig_obs.Theta,
        colormap =fig_obs.Color_scaled)
    end

    
    if Parameters.chemotaxis_flag == true
        hm = heatmap!(ax2,Parameters.pa_ph.xc,Parameters.pa_ph.yc,fig_obs.Grid_gradient,colormap = (:viridis,0.6))
    end 

    """
    # colorbar on right side of plot for the angle of the arrows
    Colorbar(f[2, 2], limits = (-1*π, 1*π), colormap = Plot_parameter.Color_arrows,
        flipaxis = true,label = "angle [rad]")
    """
    Colorbar(f[2, end+1], hm, label = "concentration [ ]")

    return ax1,ax2,fig_obs
end



function animation_step!(ax1,ax2,agent_list::CircularBuffer{Vector{agent}},Grids::Tuple{Matrix{Int64}, 
    Matrix{Float64}},fig_obs::Figure_observables, Plot_parameter::plot_parameters,Parameters::parameters)
    
    fig_obs.Time[] = fig_obs.Time[]+1
    fig_obs.Time_physical[] = fig_obs.Time_physical[] + Parameters.pa_ph.Δt_walker

    grid = Grids[1]
    fig_obs.Grid_contour[] = replace(u -> (u>0) ? 0 : u, grid)
    
    # only update scalar type of  core func if toogle is passive 
    if Plot_parameter.type_cor_func == "scalar"
        
        push!(fig_obs.V_a[], Point2f(fig_obs.Time_physical[],(Plot_parameter.cor_func_scalar(agent_list[end],Grids,Parameters))))
        
        fig_obs.V_a[] = fig_obs.V_a[]
        if length(fig_obs.V_a[]) > 3
            autolimits!(ax1)
            xlims!(ax1, fig_obs.V_a[][1][1],fig_obs.V_a[][end][1])
        end 

    elseif Plot_parameter.type_cor_func == "vector" 
        dist_vec =  Plot_parameter.cor_func_vector(agent_list[end],Grids, Parameters)
        hans = CircularBuffer{Point2f}(length(dist_vec))
        append!(hans, Point2f.(1:length(dist_vec),dist_vec))
        fig_obs.V_a[] =  hans
        
        autolimits!(ax1)
    end


    fig_obs.Grid_gradient[] = Grids[2]

    fig_obs.X.val = [i.x_pos/Parameters.arrow_to_grid_fac_x for i in agent_list[end]]
    fig_obs.Y.val = [i.y_pos/Parameters.arrow_to_grid_fac_y for i in agent_list[end]]
    fig_obs.U.val = [i.dir_x for i in agent_list[end]]
    fig_obs.V.val = [i.dir_y for i in agent_list[end]]
    Makie.update!(ax2.scene.plots[2],fig_obs.Y[], fig_obs.X[], fig_obs.V[], fig_obs.U[], color = fig_obs.Theta[] )
    fig_obs.V[] = [i.dir_y for i in agent_list[end]]
    
end


function animation_step!(ax1,agent_list::Vector{Vector{agent}},Grid_vec_sparse::Vector{SparseArrays.SparseMatrixCSC{Int64, Int64}},
    grids::Tuple{Matrix{Int64}, Matrix{Float64}},fig_obs::Figure_observables,
    Plot_parameter::plot_parameters, Parameters::parameters)
   
    fig_obs.Time[] = fig_obs.Time[]+1
    fig_obs.Time_physical[] = fig_obs.Time_physical[] + Parameters.pa_ph.Δt_walker
    
    grid = grids[1]
    reconstruct_grid!(Grid_vec_sparse,fig_obs.Time[],grid)
    fig_obs.Grid_contour[] = grid

    # ax limit cannot be controlled by observables, has to be ajusted manually in every update step
    # only update scalar type of  core func if toogle is passive 
    if Plot_parameter.type_cor_func == "scalar"
      
        
        push!(fig_obs.V_a[], Point2f(fig_obs.Time_physical[],(Plot_parameter.cor_func_scalar(agent_list[fig_obs.Time[]],grids,Parameters))))
        fig_obs.V_a[] = fig_obs.V_a[]
        if length(fig_obs.V_a[]) > 3
            autolimits!(ax1)
            xlims!(ax1, fig_obs.V_a[][1][1],fig_obs.V_a[][end][1])
        end 
    
    elseif Plot_parameter.type_cor_func == "vector" 
        
        dist_vec =  Plot_parameter.cor_func_vector(agent_list[fig_obs.Time[]],grids, Parameters)
        hans = CircularBuffer{Point2f}(length(dist_vec))
        append!(hans, Point2f.(1:length(dist_vec),dist_vec))
        fig_obs.V_a[] =  hans
        autolimits!(ax1)
    end

    #update vector plot , arrow color gets updated implicitly as it depends on u and v ;) 
    fig_obs.X.val = [i.x_pos/Parameters.arrow_to_grid_fac_x for i in agent_list[fig_obs.Time[]]]
    fig_obs.Y.val = [i.y_pos/Parameters.arrow_to_grid_fac_y for i in agent_list[fig_obs.Time[]]]
    fig_obs.U.val = [i.dir_x for i in agent_list[fig_obs.Time[]]]
    fig_obs.V[] = [i.dir_y for i in agent_list[fig_obs.Time[]]]
    
end


"""
    res_scaling(image_number; factor::Real = 3.0, plots::Int = 1)

Calculate the scaled width and height of an image based on the given factor and number of plots.

# Arguments
- `image_number::Int`: The number of the image to be scaled.
- `factor::Real`: The scaling factor for the image dimensions. Default is 3.0.
- `plots::Int`: The number of plots to be considered in the width calculation. Default is 1.

# Returns
- `width::Int`: The scaled width of the image.
- `height::Int`: The scaled height of the image.

# Example
```julia
width, height = res_scaling(10, factor=2.5, plots=2)
```
"""
function res_scaling(image_number; factor::Real = 3.0, plots::Int = 1, plots_per_row::Int = 5)
    # Initialize the counter
    c = image_number *plots

    # Calculate the scaled width and height of the image
    width = round(Int64, factor * 1000) 
    height = round(Int64, factor * 1000/plots_per_row * (c ÷ plots_per_row + ((c % plots_per_row == 0) ? 0 : 1)))

    return width, height
end