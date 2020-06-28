mutable struct TmdForcing
 
    X_wr   :: AbstractArray{Float64, 4}
    T_wr   :: AbstractArray{Float64, 3}
    S_wr   :: AbstractArray{Float64, 3}
 
    XSAS_wr   :: AbstractArray{Float64, 3}
    TSAS_wr   :: AbstractArray{Float64, 2}
    SSAS_wr   :: AbstractArray{Float64, 2}

    u_T    :: AbstractArray{Float64, 3}
    v_T    :: AbstractArray{Float64, 3}

    u_U    :: AbstractArray{Float64, 3}
    v_V    :: AbstractArray{Float64, 3}
    w_W    :: AbstractArray{Float64, 3}

    qflx2atm     :: AbstractArray{Float64, 2}
    qflx2atm_pos :: AbstractArray{Float64, 2}
    qflx2atm_neg :: AbstractArray{Float64, 2}

    qflx_X :: AbstractArray{Float64, 3}
    qflx_T :: AbstractArray{Float64, 2}
    qflx_S :: AbstractArray{Float64, 2}

    τx       :: AbstractArray{Float64, 2}
    τy       :: AbstractArray{Float64, 2}

    nswflx :: AbstractArray{Float64, 2}
    swflx  :: AbstractArray{Float64, 2}
    vsflx  :: AbstractArray{Float64, 2}
    ifrac  :: AbstractArray{Float64, 2}

    h_ML   :: AbstractArray{Float64, 2}


    function TmdForcing(
        env :: TmdEnv,
    )
        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz
        NX = env.NX

        X_wr = allocate(Float64, Nz, Nx, Ny, NX)
        T_wr = view(X_wr, :, :, :, 1)
        S_wr = view(X_wr, :, :, :, 2)
 
        XSAS_wr = allocate(Float64, Nx, Ny, NX)
        TSAS_wr = view(XSAS_wr, :, :, 1)
        SSAS_wr = view(XSAS_wr, :, :, 2)
 
        u_T = allocate(Float64, Nz, Nx, Ny)
        v_T = allocate(Float64, Nz, Nx, Ny)
 
        u_U = allocate(Float64, Nz, Nx, Ny)
        v_V = allocate(Float64, Nz, Nx, Ny+1)
        w_W = allocate(Float64, Nz+1, Nx, Ny)
 
        qflx2atm = allocate(Float64, Nx, Ny)
        qflx2atm_pos = allocate(Float64, Nx, Ny)
        qflx2atm_neg = allocate(Float64, Nx, Ny)

        qflx_X  = allocate(Float64, Nx, Ny, NX)
        qflx_T = view(qflx_X, :, :, 1)
        qflx_S = view(qflx_X, :, :, 2)
 
        τx = allocate(Float64, Nx, Ny)
        τy = allocate(Float64, Nx, Ny)
        
        nswflx = allocate(Float64, Nx, Ny)
        swflx  = allocate(Float64, Nx, Ny)
        vsflx  = allocate(Float64, Nx, Ny)
        ifrac  = allocate(Float64, Nx, Ny)
        
        h_ML = allocate(Float64, Nx, Ny)
 
         return new(
            X_wr,
            T_wr,
            S_wr,
         
            XSAS_wr,
            TSAS_wr,
            SSAS_wr,

            u_T,
            v_T,

            u_U,
            v_V,
            w_W,

            qflx2atm,
            qflx2atm_pos,
            qflx2atm_neg,

            qflx_X,
            qflx_T,
            qflx_S,

            τx,
            τy,

            nswflx,
            swflx,
            vsflx,
            ifrac,

            h_ML,
        )
    end
end


