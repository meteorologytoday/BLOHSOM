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


    function DynState(env)

        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz

        u_total =  zeros(Float64, Nx, Ny, Nz)
        v_total =  zeros(Float64, Nx, Ny, Nz)

        Φ =  zeros(Float64, Nx, Ny)
        
        U =  zeros(Float64, Nx, Ny  )
        V =  zeros(Float64, Nx, Ny)


        u =  zeros(Float64, Nx, Ny, Nz)
        v =  zeros(Float64, Nx, Ny, Nz)
        b =  zeros(Float64, Nx, Ny, Nz)

        G_u = zeros(Float64, Nx, Ny, Nz, 3)
        G_v = zeros(Float64, Nx, Ny, Nz, 3)

        return new(
            u_total, v_total,
            Φ, U, V,
               u, v,
            G_u, G_v,
        )
        
    end
end



