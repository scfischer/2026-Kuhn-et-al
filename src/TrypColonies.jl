module TrypColonies

using Distributions, LinearAlgebra, SparseArrays
import Random as r
import Serialization as s
using DataStructures, Parameters
using DataFrames, ColorSchemes, GLMakie, ImageFiltering
#using Polyester
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
using FFTW
    
export  @h,
        create_grids,
        parameters,
        plot_parameters,
        order_parameter,
        initialize_system,
        make_sys_interactive!,
        makie_config!,
        ini_animation,
        create_buttons,
        update_directions!,
        update_position!,
        strenghten_boundary!,
        animation_step!,
        update_grid!,
        make_para_dict_list,
        sweep_and_save_parallel_t,
        make_para_dict_list,
        create_para_struct,
        order_parameter_circle,
        reconstruct_grid_from_scratch,
        find_data_path,
        foldersize,
        create_data_path_vector,
        create_para_struct_vec,
        ana_order_para_vicsek,
        create_sweep_para_vec,
        ana_order_para_circle,
        create_exp_scalar_grid,
        angle_between_vec,
        divison!,
        parameters_physical,
        create_diff_grids,
        diffusion_2D!,
        adsorb!,
        radial_distribution,
        angular_metric,
        pair_cor_metric, 
        order_parameter_perpendic,
        remove_fields_dict!,
        sweep_and_save_parallel_t_static,
        initialize_system_static,
        b_w_im,
        b_w_oc,
        res_scaling,
        area_gain,
        transform_data,
        radial_density_c,
        adsorbed_material,
        cell_gain,
        og_size,
        fourier_ik,
        roughness,
        scale_color_map




include("Agent_grid.jl")
include("Position.jl")
include("Direction.jl")
include("Analysis.jl")
include("Sweep_functions.jl")
include("Visu.jl")
include("Auxiliary.jl")
include("Division.jl")
include("Diffusion.jl")
include("Adsorption.jl")
    
    

end
