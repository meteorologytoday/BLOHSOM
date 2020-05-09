mutable struct DynForcing
    
    B        :: AbstractArray{Float64, 3}

    # On T grid, wind stress in real east (x) and north (y)
    τx_raw  :: AbstractArray{Float64, 2} 
    τy_raw  :: AbstractArray{Float64, 2}  

    # wind stress along ocean x, y coordinate
    τx  :: AbstractArray{Float64, 2}
    τy  :: AbstractArray{Float64, 2}

    function DynForcing(env)

        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz

        B       =  zeros(Float64, Nx, Ny, Nz)

        τx_raw = zeros(Float64, Nx, Ny)
        τy_raw = zeros(Float64, Nx, Ny)

        τx = zeros(Float64, Nx, Ny)
        τy = zeros(Float64, Nx, Ny+1)



        return new(
            B,
            τx_raw,
            τy_raw,
            τx,
            τy,
        )
        
    end
end



