mutable struct DynState
    
    # total velocity, buoyancy
    u_total  :: AbstractArray{Float64, 3}
    v_total  :: AbstractArray{Float64, 3}

    Φ  :: AbstractArray{Float64, 2}
    
    # barotropic component
    U  :: AbstractArray{Float64, 2}
    V  :: AbstractArray{Float64, 2}

    # Baroclinic component
    u  :: AbstractArray{Float64, 3}
    v  :: AbstractArray{Float64, 3}

    # Excpetion 2
    G_u :: AbstractArray{Float64, 4}  # G term for baroclinic mode.   (Nc_c, Nx+1, Ny, time)
    G_v :: AbstractArray{Float64, 4}  # G term for baroclinic mode.   (Nc_c, Nx, Ny+1, time)


    # Forcing variables #
    B        :: AbstractArray{Float64, 3}

    # On T grid, wind stress in real east (x) and north (y)
    τx_raw  :: AbstractArray{Float64, 2} 
    τy_raw  :: AbstractArray{Float64, 2}  

    # wind stress along ocean x, y coordinate
    τx  :: AbstractArray{Float64, 2}
    τy  :: AbstractArray{Float64, 2}

    function DynState(env)

        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz

        u_total =  allocate(Float64, Nx, Ny, Nz)
        v_total =  allocate(Float64, Nx, Ny+1, Nz)

        Φ =  allocate(Float64, Nx, Ny)
        
        U =  allocate(Float64, Nx, Ny  )
        V =  allocate(Float64, Nx, Ny+1)


        u =  allocate(Float64, Nx, Ny, Nz)
        v =  allocate(Float64, Nx, Ny+1, Nz)
        b =  allocate(Float64, Nx, Ny, Nz)

        G_u = allocate(Float64, Nx, Ny, Nz, 3)
        G_v = allocate(Float64, Nx, Ny+1, Nz, 3)



        # Forcing variables

        B       =  allocate(Float64, Nx, Ny, Nz)

        τx_raw = allocate(Float64, Nx, Ny)
        τy_raw = allocate(Float64, Nx, Ny)

        τx = allocate(Float64, Nx, Ny)
        τy = allocate(Float64, Nx, Ny+1)

        return new(
            u_total, v_total,
            Φ, U, V,
               u, v,
            G_u, G_v,

            B,
            τx_raw,
            τy_raw,
            τx,
            τy,
        )
        
    end
end



