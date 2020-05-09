function decompose!(;
    total   :: AbstractArray{Float64, 3},
    mean    :: AbstractArray{Float64, 2},
    anomaly :: AbstractArray{Float64, 3},
    va       :: VerticalAverager,
)

    Nx, Ny, Nz = size(anomaly)

    
    calAverage_c2s!(
        va,
        total,
        mean;
    )

    for i=1:Nx, j=1:Ny

        m = mean[i, j]
        for k=1:Nz
            anomaly[i, j, k] = total[i, j, k] - m
        end

    end
end

function advectDynamic!(
    model   :: DynModel,
)
    env   = model.env
    state = model.state
    core  = model.core

#=
  0.041038 seconds (106 allocations: 9.164 MiB, 2.54% gc time)
  0.045352 seconds (2.16 M allocations: 32.969 MiB, 2.54% gc time)
  0.026442 seconds (1.44 M allocations: 21.995 MiB, 5.05% gc time)
  0.005663 seconds (21 allocations: 1.572 MiB)
  0.001000 seconds (100 allocations: 3.063 KiB)
=#
    decomposeModes!(model)
#    println("Do Diffusion!")
#    println("calAuxV!")
    calAuxV!(model)

#    println("calAuxΦ!")
    calAuxΦ!(model)

#    println("solveΦ!")
    solveΦ!(model)

#    model.state.Φ .= model.env.Φ_total

#    println("updateV!")
    updateV!(model)

    decomposeModes!(model)
    doHDiffusionBarotropic!(model)
    #println(format("(u, v) = ({:.2f}, {:.2f})", state.u_total[1, 5, 5], state.v_total[1,5,5]))
end

function doHDiffusionBarotropic!(
    model :: DynModel,
)
 
    state = model.state
    core = model.core
    env = model.env
    wksp = core.wksp
    layers = core.layers 
    
    wksp_sU = getSpace!(wksp, :sU)
    wksp_sV = getSpace!(wksp, :sV)

    wksp_sU_sol = getSpace!(wksp, :sU)
    wksp_sV_sol = getSpace!(wksp, :sV)

#    mul2!(wksp_sU, core.s_ops.U_interp_T, state.U)
#    mul2!(wksp_sV, core.s_ops.V_interp_T, state.V)

#    println(core.diffusion_solver_Kh_barotropic.K)

    wksp_sU .= state.U 
    wksp_sV .= state.V

    solveDiffusion!(
        core.diffusion_solver_Kh_barotropic, :U,
        wksp_sU,
        wksp_sU_sol,
    )

    solveDiffusion!(
        core.diffusion_solver_Kh_barotropic, :V,
        wksp_sV,
        wksp_sV_sol,
    )


    # ===== [ BEGIN coastal friction on barotropic mode ] =====

    ϵΔtp1 = 1 + env.Δt / env.τ_barotropic_coastfric
    r = 1.0 / (ϵΔtp1+1.0) - 1.0

    mul2!(wksp_sV, core.c_ops.coastmask_V, wksp_sV_sol)
    @. wksp_sV_sol = wksp_sV_sol + r * wksp_sV *2 

    mul2!(wksp_sU, core.c_ops.coastmask_U, wksp_sU_sol)
    @. wksp_sU_sol = wksp_sU_sol + r * wksp_sU*2

    state.U .= wksp_sU_sol
    state.V .= wksp_sV_sol

    # ===== [ END coastal friction on barotropic mode ] =====

    # baroclinic

    for k=1:env.Nz

        wksp_sU .= layers.u[k] 
        wksp_sV .= layers.v[k] 

        solveDiffusion!(
            core.diffusion_solver_Kh_baroclinic, :U,
            wksp_sU,
            layers.u[k],
        )

        solveDiffusion!(
            core.diffusion_solver_Kh_barotropic, :V,
            wksp_sV,
            layers.v[k],
        )

    end

    for k = 1:env.Nz
        @. layers.u_total[k] = state.U + layers.u[k]
        @. layers.v_total[k] = state.V + layers.v[k]
    end
end

function decomposeModes!(
    model :: DynModel
)
 
    state = model.state
    core = model.core
    c_ops = core.c_ops
    s_ops = core.s_ops
    va = core.va
    env = model.env
    wksp = core.wksp
  
    
    # cal barotropic and baroclinic components
    decompose!(total=state.u_total, mean=state.U, anomaly=state.u, va=va)
    decompose!(total=state.v_total, mean=state.V, anomaly=state.v, va=va)

end


function calAuxV!(
    model   :: DynModel,
)

    @unpack model

    state = model.state
    core = model.core
    c_ops = core.c_ops
    s_ops = core.s_ops
    va = core.va
    env = model.env
    wksp = core.wksp
  
  
    # cal τx_acc, τy_acc
    τx_acc = getSpace!(wksp, :sU)
    τy_acc = getSpace!(wksp, :sV)

    

    @. τx_acc = fr.τx / (ρ_fw * ev.H[1])
    @. τy_acc = fr.τy / (ρ_fw * ev.H[1])
 
    # cal ∇b
    ∂B∂x   = core.∂B∂x
    ∂B∂y   = core.∂B∂y

#    BBB = getSpace!(wksp, :U)
#    mul3!(BBB, c_ops.U_interp_T, fr.B)

    #println(fr.B[21:30, 45, 2])
    #println(BBB[21:30, 45, 2])
    #readline()

#    mul3!(core.∂B∂x, c_ops.T_∂x_T, fr.B)
#    mul3!(core.∂B∂y, c_ops.T_∂y_T, fr.B)

    mul3!(core.∂B∂x, c_ops.U_∂x_T, fr.B)
    mul3!(core.∂B∂y, c_ops.V_∂y_T, fr.B)


    # cal Coriolis force
    fu   = getSpace!(wksp, :V)
    fv   = getSpace!(wksp, :U)
    
#    println("size of fu: ", size(fu))
#    println("size of V_f_U: ", size(c_ops.V_f_U))
#    println("sizeof u_total: ", size(state.u_total))

    mul3!(fu, c_ops.V_f_U, state.u_total)   # fu on V grid
    mul3!(fv, c_ops.U_f_V, state.v_total)   # fv on U grid
#    mul3!(fu, c_ops.V_f_U, state.u)   # fu on V grid
#    mul3!(fv, c_ops.U_f_V, state.v)   # fv on U grid

    #println("fu: ", fu[30, 29, 1])   

 
    # ===== [ BEGIN cal (v⋅∇)v ] =====
    #=
    ∂u∂x = getSpace!(wksp, :U)
    ∂u∂y = getSpace!(wksp, :U)
    ∂v∂x = getSpace!(wksp, :V)
    ∂v∂y = getSpace!(wksp, :V)
    mul3!(∂u∂x, c_ops.U_∂x_U, state.u)
    mul3!(∂u∂y, c_ops.U_∂y_U, state.u)
    mul3!(∂v∂x, c_ops.V_∂x_V, state.v)
    mul3!(∂v∂y, c_ops.V_∂y_V, state.v)
    
    # interpolate first then multiply by ∇v
    
    # On U grid

    
    u∂u∂x = getSpace!(wksp, :U)
    v∂u∂y = getSpace!(wksp, :U)
    mul3!(v∂u∂y, c_ops.U_interp_V, state.v_total)  # store interpolated v into v∂u∂y
    for i=1:env.Nx, j=1:env.Ny, k=1:env.Nz
        u∂u∂x[i, j, k]  = state.u[i, j, k] * ∂u∂x[i, j, k]
        v∂u∂y[i, j, k] *=                    ∂u∂y[i, j, k]
    end

    # On V grid
    u∂v∂x = getSpace!(wksp, :V)
    v∂v∂y = getSpace!(wksp, :V)
    mul3!(u∂v∂x, c_ops.V_interp_U, state.u_total)  # store interpolated u into u∂v∂x
    for i=1:env.Nx, j=1:env.Ny+1, k=1:env.Nz
        u∂v∂x[i, j, k] *=                     ∂v∂x[i, j, k]
        v∂v∂y[i, j, k]  = state.v[i, j, k]  * ∂v∂y[i, j, k]
    end
    =#
    # ===== [ END cal (v⋅∇)v ] =====

    # cal G
    G_idx = core.G_idx
    Δt0 = G_idx[:now]
    Δt1 = G_idx[:one_Δt_ago]
    Δt2 = G_idx[:two_Δt_ago]


    G_u = view(state.G_u, :, :, :, Δt0)
    G_v = view(state.G_v, :, :, :, Δt0)

    @. G_u = ∂B∂x + fv
#    G_u .-= u∂u∂x
#    G_u .-= v∂u∂y
 
    @. G_v = ∂B∂y - fu
#    G_v .-= u∂v∂x
#    G_v .-= v∂v∂y


    # surface
#    G_u[:, :, 1] .+= τx_acc
#    G_v[:, :, 1] .+= τy_acc
    #println("G_u")

    #=
    println("U")
    println(state.U)

    println("u first layer")
    println(state.u[1, :, :])

    println("dudx")
    println(∂u∂x[1, :, :])
    =#

    # calculate auxiliary velocity
    Δt = env.Δt 
#    G_u = state.G_u
#    G_v = state.G_v

    G_u_Δt0 = view(state.G_u, :, :, :, Δt0)
    G_u_Δt1 = view(state.G_u, :, :, :, Δt1)
    G_u_Δt2 = view(state.G_u, :, :, :, Δt2)

    G_v_Δt0 = view(state.G_v, :, :, :, Δt0)
    G_v_Δt1 = view(state.G_v, :, :, :, Δt1)
    G_v_Δt2 = view(state.G_v, :, :, :, Δt2)


    @. core.u_aux = state.u_total + Δt * ABIII(G_u_Δt0, G_u_Δt1, G_u_Δt2)
    @. core.v_aux = state.v_total + Δt * ABIII(G_v_Δt0, G_v_Δt1, G_v_Δt2)

    #=
    @time for i=1:env.Nx, j=1:env.Ny, k=1:env.Nz
        core.u_aux[i, j, k] = state.u_total[i, j, k] + Δt *
           ABIII(
                state.G_u[i, j, k, Δt0],
                state.G_u[i, j, k, Δt1],
                state.G_u[i, j, k, Δt2],
           )
    end

    #println("G_u")
    #println(state.G_u[1, 2, 2, Δt0])

    for i=1:env.Nx, j=1:env.Ny+1, k=1:env.Nz
        core.v_aux[i, j, k] = state.v_total[i, j, k] + Δt *
            ABIII(
                G_v[i, j, k, Δt0],
                G_v[i, j, k, Δt1],
                G_v[i, j, k, Δt2],
            )
    end
    =#

    G_idx[:now] = Δt2
    G_idx[:one_Δt_ago] = Δt0
    G_idx[:two_Δt_ago] = Δt1

end

function calAuxΦ!(
    model :: DynModel,
)
    # cal mean aux_v
    core = model.core
    env  = model.env
    state = model.state

    s_ops = core.s_ops
    va    = core.va
    wksp  = core.wksp
    mask  = env.mask
    Δt    = env.Δt

 
#    u_T     = getSpace!(wksp, :sT)
#    v_T     = getSpace!(wksp, :sT)
   
    Φu     = getSpace!(wksp, :sU)
    Φv     = getSpace!(wksp, :sV)
    DIV_Φu = getSpace!(wksp, :sT)
    DIV_Φv = getSpace!(wksp, :sT)

    calAverage_c2s!(va, core.u_aux, Φu)
    calAverage_c2s!(va, core.v_aux, Φv)

#    mul2!(Φu, s_ops.U_interp_T, u_T)
#    mul2!(Φv, s_ops.V_interp_T, v_T)

    @. Φu *= env.Φ_total
    @. Φv *= env.Φ_total

    mul2!(DIV_Φu, s_ops.T_DIVx_U,   Φu)
    mul2!(DIV_Φv, s_ops.T_DIVy_V,   Φv)
 
    @. core.Φ_aux = state.Φ - Δt * (DIV_Φu + DIV_Φv)

end

function solveΦ!(
    model :: DynModel,
)

    env    = model.env
    core   = model.core
    state  = model.state
    solver = core.Φ_solver
    
    #solvePhi!(core.Φ_solver, core.Φ_aux, state.Φ)
   
    #return

    lhs = getSpace!(core.wksp, :sT)
    rhs = getSpace!(core.wksp, :sT)

    lhs_eT = view(getSpace!(core.wksp, :sT), 1:solver.eT_length)
    rhs_eT = view(getSpace!(core.wksp, :sT), 1:solver.eT_length)

    α = solver.α
    
    @. rhs = - core.Φ_aux * α

    mul!(rhs_eT, solver.eT_send_T, view(rhs, :))
 
    ldiv!(
        lhs_eT,
        solver.tool_mtx.MoLap,
        rhs_eT,
    )
        
    mul!(view(state.Φ, :), solver.T_send_eT, lhs_eT)
end

function updateV!(
    model :: DynModel,
)
   
    core  = model.core
    env   = model.env
    state = model.state

    s_ops = core.s_ops
    wksp = core.wksp
    Δt   = model.env.Δt
    
    Δt∂Φ∂x = getSpace!(wksp, :sU)
    Δt∂Φ∂y = getSpace!(wksp, :sV)
     
    mul2!(Δt∂Φ∂x, s_ops.U_∂x_T, state.Φ)
    mul2!(Δt∂Φ∂y, s_ops.V_∂y_T, state.Φ)

    Δt∂Φ∂x .*= Δt
    Δt∂Φ∂y .*= Δt

    # TODO: NEED to use AUXILIRAY U, V instead of old U, V 
    ΔtfricU = getSpace!(wksp, :sU)
    ΔtfricV = getSpace!(wksp, :sV)

    r = Δt/env.τ_barotropic_bottomfric
    @. ΔtfricU = state.U * r
    @. ΔtfricV = state.V * r

    for k=1:env.Nz
        u_total = view(state.u_total, :, :, k)
        v_total = view(state.v_total, :, :, k)
        u_aux = view(core.u_aux, :, :, k)
        v_aux = view(core.v_aux, :, :, k)

        @. u_total = u_aux - Δt∂Φ∂x - ΔtfricU
        @. v_total = v_aux - Δt∂Φ∂y - ΔtfricV
    end

#= 
    for i=1:env.Nx, j=1:env.Ny, k=1:env.Nz
        state.u_total[i, j, k] = core.u_aux[i, j, k] - Δt∂Φ∂x[i, j] - ΔtfricU[i, j]
        state.v_total[i, j, k] = core.v_aux[i, j, k] - Δt∂Φ∂y[i, j] - ΔtfricV[i, j]
    end
 =#
    #projVertical_c2f!(core.va, state.u_total, state.u_f)
    #projVertical_c2f!(core.va, state.v_total, state.v_f)
   
end

