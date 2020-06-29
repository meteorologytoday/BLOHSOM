function determineVelocity!(
    m :: TmdModel,
)
    @fast_extract m

    wksp  = co.wksp

    div   = getSpace!(wksp, :T)
    tmp_T = getSpace!(wksp, :T)

#    u_U = getSpace!(wksp, :U)
#    v_V = getSpace!(wksp, :V)


#    u_U .= st.u_U
#    v_V .= st.v_V

#    mul_autoflat!(u_U, co.ASUM.filter_U, st.u_U)
#    mul_autoflat!(v_V, co.ASUM.filter_V, st.v_V)
 
    calDIV!(
        ASUM = co.ASUM,
        u_U = st.u_U,
        v_V = st.v_V,
        div   = div,
        workspace = tmp_T
    )

    calVerVelBnd!(
        gi    = ev.gi,
        Nx    = ev.Nx,
        yrng  = ev.update_yrng_T,
        Nz    = ev.Nz_av,
        w_bnd = st.w_W,
        hs    = co.Δz_T,
        div   = div,
        mask3 = ev.mask3,
    )

end

# This function is deriveing the local 
# tendency without updating the state
# so that parallization is possible.
function advectTracer_part1!( 
    m   :: TmdModel,
)

    @fast_extract m

    calFLDOPartition!(m)

#    println("update_yrng_T: ", ev.update_yrng_T)

    for x = 1:ev.NX

        ΔX   = view(st.ΔX,      :, :, x)
        X_ML = view(st.X_ML,    :, :, x)
        X    = view(co.cols.X,  :, :, x)


        for i=1:ev.Nx, j=ev.update_yrng_T
            ΔX[i, j] = mixFLDO!(
                qs   = X[i, j],
                zs   = co.cols.z_bnd_av[i, j],
                hs   = co.cols.Δz_T[i, j],
                q_ML = X_ML[i, j],
                FLDO = st.FLDO[i, j],
                FLDO_ratio_top = st.FLDO_ratio_top[i, j],
                FLDO_ratio_bot = st.FLDO_ratio_bot[i, j],
            )
        end

    end

#    println("AFTER:T[1:2]", st.X[1:2, 3, 3, 1], ";h_ML=", st.h_ML[3, 3])

    # Pseudo code
    # 1. calculate tracer flux
    # 2. calculate tracer flux divergence

    #T_sum_old = sum(co.ASUM.T_Δvol_T * view(st.T, :))

    calDiffAdv_QUICKEST_SpeedUp!(m, ev.Δt_substep)

end
 
function advectTracer_part2!( 
    m   :: TmdModel,
)

    @fast_extract m

#=
    # deal with top layer
    w_sfc    = view(st.w_W,  1, :, :)
    Δz_T_sfc = ev.z_bnd[1] - ev.z_bnd[2]#view(co.Δz_T, 1, :, :)
    for x = 1:ev.NX
        X = view(st.X, 1, :, :, x) 
#        @. view(co.XFLUX_top, :, :, x) = X * w_sfc
        @. X -= Δt * w_sfc * X / Δz_T_sfc
    end
    #st.Xsfcsponge += co.XFLUX_top * Δt 
=#
    X          = view(st.X,          :, :, ev.update_yrng_T, :)
    XFLUX_CONV = view(co.XFLUX_CONV, :, :, ev.update_yrng_T, :)
    ΔX         = view(st.ΔX,    :, ev.update_yrng_T, :)
    dΔXdt      = view(st.dΔXdt, :, ev.update_yrng_T, :)

    @. X   += ev.Δt_substep * XFLUX_CONV
    @. ΔX  += ev.Δt_substep * dΔXdt

#    println("co.XFLUX_CONV * Δt = ", Δt * co.XFLUX_CONV[:, 60, 59, 1])
#    println(" u_U = ", st.u_Uco.XFLUX_CONV[:, 60, 59, 1])

    #T_sum_new = sum(co.ASUM.T_Δvol_T * view(st.T, :))
    #println("Total_heat_new: ", T_sum_new, "; change: ", (T_sum_new-T_sum_old)/T_sum_old * 100, "%")


#    T_sum_new = sum(co.ASUM.T_Δvol_T * view(st.T, :))

#    println("Tsum: ", T_sum_old, " => ", T_sum_new, "; change = ", (T_sum_new - T_sum_old) / T_sum_old * 100.0, " %")


#    println("## Before T_ML=", st.X_ML[3,3,1])
    for x = 1:ev.NX
        ΔX   = view(st.ΔX,      :, :, x)
        X_ML = view(st.X_ML,    :, :, x)
        X    = view(co.cols.X,  :, :, x)

        for i=1:ev.Nx, j=ev.update_yrng_T
            X_ML[i, j] = unmixFLDOKeepDiff!(
                qs   = X[i, j],
                zs   = co.cols.z_bnd_av[i, j],
                hs   = co.cols.Δz_T[i, j],
                h_ML = st.h_ML[i, j],
                FLDO = st.FLDO[i, j],
                Nz   = ev.Nz_av[i, j],
                Δq   = ΔX[i, j],
            )
        end
    end

end

function calVerVelBnd!(;
    gi       :: PolelikeCoordinate.GridInfo,
    Nx       :: Integer,
    yrng     :: UnitRange,
    Nz       :: AbstractArray{Int64, 2},
    w_bnd    :: AbstractArray{Float64, 3},   # ( Nz_bone+1, Nx  , Ny   )
    hs       :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    div      :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
    mask3    :: AbstractArray{Float64, 3},   # ( Nz_bone  , Nx  , Ny   )
)

#    println(size(w_bnd))
#    println(size(div))
#    println(size(hs))
#    local tmp = tmp_σ = 0.0
    for i=1:Nx, j=yrng
        
        _Nz = Nz[i, j]
        w_bnd[1, i, j]     = 0.0
        w_bnd[_Nz+1, i, j] = 0.0

        for k=_Nz:-1:1

            if mask3[k, i, j] == 0.0
                break
            end
            
            #w_bnd[k+1, i, j] = w_bnd[k, i, j] + div[k, i, j] * hs[k, i, j]
            w_bnd[k, i, j] = w_bnd[k+1, i, j] - div[k, i, j] * hs[k, i, j]
        end
        
        #=
        for k=1:_Nz

            if mask3[k, i, j] == 0.0
                break
            end
            
            w_bnd[k+1, i, j] = w_bnd[k, i, j] + div[k, i, j] * hs[k, i, j]
        end
        =#

#        tmp   += w_bnd[Nz[i, j]+1, i, j] * gi.dσ[i, j]
#        tmp_σ += gi.dσ[i, j]
    end

#   #println("tmp: ", tmp, "; tmp_σ: ", tmp_σ, "; Average w: ", tmp/tmp_σ)

end



