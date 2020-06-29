mutable struct TmdState
    
    # Prognostic variables 
    X         :: AbstractArray{Float64, 4}
    T         :: AbstractArray{Float64, 3}
    S         :: AbstractArray{Float64, 3}
  
    X_ML     :: AbstractArray{Float64, 3}
    T_ML     :: AbstractArray{Float64, 2}
    S_ML     :: AbstractArray{Float64, 2}

    
    ΔX        :: AbstractArray{Float64, 3}
    ΔT        :: AbstractArray{Float64, 2}
    ΔS        :: AbstractArray{Float64, 2}

    dΔXdt     :: AbstractArray{Float64, 3}
    dΔTdt     :: AbstractArray{Float64, 2}
    dΔSdt     :: AbstractArray{Float64, 2}


    h_ML      :: AbstractArray{Float64, 2}
    h_MO      :: AbstractArray{Float64, 2}
    
    FLDO           :: AbstractArray{Int64,     2}
    FLDO_ratio_top :: AbstractArray{Float64,   2}
    FLDO_ratio_bot :: AbstractArray{Float64,   2}

    b        :: AbstractArray{Float64, 3}
    b_ML     :: AbstractArray{Float64, 2}
 
    b_mixed  :: AbstractArray{Float64, 3}
    B        :: AbstractArray{Float64, 3}

    # Diagnostic variables
    intX             :: AbstractArray{Float64, 3} # Total X content, vertical integrated quantity
    dintXdt          :: AbstractArray{Float64, 3}

    intT             :: AbstractArray{Float64, 2} # Total heat content
    dintTdt          :: AbstractArray{Float64, 2} # Total heat content change rate
    intS             :: AbstractArray{Float64, 2} # Total salt
    dintSdt          :: AbstractArray{Float64, 2} # Total salt change rate

    qflx_X_correction :: AbstractArray{Float64, 3}
    qflx_T_correction :: AbstractArray{Float64, 2}
    qflx_S_correction :: AbstractArray{Float64, 2}

    XFLUX_DIV_implied :: AbstractArray{Float64, 3}
    TFLUX_DIV_implied :: AbstractArray{Float64, 2}
    SFLUX_DIV_implied :: AbstractArray{Float64, 2}


    # Forcing variables
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

    pres_h_ML   :: AbstractArray{Float64, 2} # prescribed h_ML


    function TmdState()
        return new()
    end

    function TmdState(env :: TmdEnv)

        return TmdState(
            env.Nx,
            env.Ny,
            env.Nz,
            env.NX,
        )

    end
 
    function TmdState(
        Nx :: Int64,
        Ny :: Int64,
        Nz :: Int64,
        NX :: Int64,
    )

        # Will be using dimension length to
        # determine the dimension label.
        # This is only for convenience in the
        # `CreateView` method.
        diffs = [
            Ny-Nx,
            Ny-Nz,
            Ny-NX,
        ]

        for diff in diffs
            if abs(diff) < 2
                throw(ErrorException("`Ny` must be different to other dimensions by greater or equal to 2."))
            end
        end



        #

        X = allocate(Float64, Nz, Nx, Ny, NX)
        T = view(X, :, :, :, 1)
        S = view(X, :, :, :, 2)


        X_ML = allocate(Float64, Nx, Ny, NX)
        T_ML = view(X_ML, :, :, 1)
        S_ML = view(X_ML, :, :, 2)
        
        ΔX = allocate(Float64, Nx, Ny, NX)
        ΔT = view(ΔX, :, :, 1)
        ΔS = view(ΔX, :, :, 2)
 
        dΔXdt = allocate(Float64, Nx, Ny, NX)
        dΔTdt = view(dΔXdt, :, :, 1)
        dΔSdt = view(dΔXdt, :, :, 2)

       
        h_ML = allocate(Float64, Nx, Ny)
        h_MO = allocate(Float64, Nx, Ny)
        
        FLDO = allocate(Int64, Nx, Ny)
        FLDO_ratio_top = allocate(Float64, Nx, Ny)
        FLDO_ratio_bot = allocate(Float64, Nx, Ny)

        b    = allocate(Float64, Nz, Nx, Ny)
        b_ML = allocate(Float64, Nx, Ny)
        
        b_mixed = allocate(Float64, Nz, Nx, Ny)
        B       = allocate(Float64, Nz, Nx, Ny)

        #

        intX = allocate(Float64, Nx, Ny, NX)
        intT = view(intX, :, :, 1)
        intS = view(intX, :, :, 2)

        dintXdt = allocate(Float64, Nx, Ny, NX)
        dintTdt = view(dintXdt, :, :, 1)
        dintSdt = view(dintXdt, :, :, 2)
 
        qflx_X_correction = allocate(Float64, Nx, Ny, NX)
        qflx_T_correction = view(qflx_X_correction, :, :, 1)
        qflx_S_correction = view(qflx_X_correction, :, :, 2)

        XFLUX_DIV_implied = allocate(Float64, Nx, Ny, NX)
        TFLUX_DIV_implied = view(XFLUX_DIV_implied, :, :, 1)
        SFLUX_DIV_implied = view(XFLUX_DIV_implied, :, :, 2)


        #

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
        
        pres_h_ML = allocate(Float64, Nx, Ny)


        return new(
            X,
            T,
            S,
           
            X_ML,
            T_ML,
            S_ML,
           
            ΔX,
            ΔT,
            ΔS,

            dΔXdt,
            dΔTdt,
            dΔSdt,

            h_ML,
            h_MO,
            
            FLDO,
            FLDO_ratio_top,
            FLDO_ratio_bot,

            b,
            b_ML,

            b_mixed,
            B,

            intX,
            dintXdt,

            intT,
            dintTdt,
            intS,
            dintSdt,

            qflx_X_correction,
            qflx_T_correction,
            qflx_S_correction,

            XFLUX_DIV_implied,
            TFLUX_DIV_implied,
            SFLUX_DIV_implied,


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

            pres_h_ML,

        )
    end

end

function createMirror(
    ref_state    :: TmdState,
    ref_Ny       :: Int64,
    yrng_T       :: UnitRange,
    yrng_V       :: UnitRange,
)

    mirror_state = TmdState()

    for varsym in fieldnames(TmdState)

        var = Core.getfield(ref_state, varsym)
        
        dims = size(var)

        _rng = []
        for (i, N) in enumerate(dims)

            if N == ref_Ny
                push!(_rng, yrng_T)
            elseif N == ref_Ny+1
                push!(_rng, yrng_V)
            else
                push!(_rng, Colon())
            end
            
        end

        Core.setfield!(mirror_state, varsym, view(var, _rng...))

    end

    return mirror_state
end
