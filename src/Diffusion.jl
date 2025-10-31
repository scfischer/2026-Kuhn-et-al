@init_parallel_stencil(Threads, Float64, 2)


@parallel function compute_flux!(qTx, qTy, T, Lam, dx, dy)
    @all(qTx) = -Lam*@d_xi(T)/dx
    @all(qTy) = -Lam*@d_yi(T)/dy
    return
end


@parallel function compute_update!(T, qTx, qTy, dt, dx, dy, decay_rate)
    @inn(T) = @inn(T) - dt*(@d_xa(qTx)/dx + @d_ya(qTy)/dy) - dt*decay_rate*@inn(T)
    return
end

function create_diff_grids(para::parameters)
    para_ph = para.pa_ph    
    qTx      = @zeros(para_ph.N[1]-1,para_ph.N[2]-2)
    qTy      = @zeros(para_ph.N[1]-2,para_ph.N[2]-1)
    return (qTx, qTy)
end

function diffusion_2D!(grids::Tuple{Matrix{Int64}, Matrix{Float64}},
    diff_grids::Tuple{Matrix{Float64}, Matrix{Float64}},para::parameters)
    
    T = grids[2]
    qTx = diff_grids[1]
    qTy = diff_grids[2]

    para_ph = para.pa_ph
    dt = para_ph.Δt_diff
    Lam  = para_ph.Diff_coe
    decay_rate = para_ph.decay_rate

    @parallel compute_flux!(qTx, qTy, T, Lam, para_ph.dx, para_ph.dy)
    @parallel compute_update!(T, qTx, qTy, dt, para_ph.dx, para_ph.dy, decay_rate)
end

