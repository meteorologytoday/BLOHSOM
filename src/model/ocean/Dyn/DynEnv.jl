mutable struct DynEnv
    
    #
    # This shallow water model is fixed height
    # and advect high vertical resultion while
    # horizontal is split into coarser grids
    #
    # c : coarse vertical grid
    # f : fine   vertical grid
    #
    # ----------------+---+---+
    # variable        | c | f |
    # ----------------+---+---+
    # u, v            | v | v |
    # b               | v |   |
    # T, S, X         |   | v |
    # H               | v |   |
    #
    # Variables are arranged in Arakawa-C grid
    # 
        
    gi :: PolelikeCoordinate.GridInfo

    Δt :: Float64
    Kh_barotropic :: Float64
    Kh_baroclinic :: Float64

    Nx :: Int64
    Ny :: Int64
    Nz :: Int64
    
    z_bnd :: AbstractArray{Float64, 1}

    H  :: AbstractArray{Float64, 1}
    H_total :: Float64
    Φ_total :: Float64
   
    mask :: AbstractArray{Float64, 2}     # where mixed layer model is active

    mode :: Symbol  # `PROG`, `EKMAN`
    # If any mode `EKMAN` takes only three layers:
    # Ekman layer, return flow layer, static layer

    ϵ :: Float64
    s :: AbstractArray{Complex{Float64}, 2}

    function DynEnv(;
        gi                  :: PolelikeCoordinate.GridInfo,
        Δt                  :: Float64,
        Kh_barotropic       :: Float64,
        Kh_baroclinic       :: Float64,
        Nx                  :: Int64,
        Ny                  :: Int64,
        z_bnd               :: AbstractArray{Float64, 1},
        mode                :: Symbol,
        mask                :: Union{AbstractArray{Float64, 2}, Nothing},
        ϵ                   :: Float64 = 0.0,
    )
 
        z_bnd = copy(z_bnd)
        Nz = length(z_bnd) - 1

        H = z_bnd[1:end-1] - z_bnd[2:end]

        # If mask is not provided
        if mask == nothing
            mask = ones(Nx, Ny)
        end

        H_total = sum(H)
        Φ_total = g * H_total
        
        println("H_total: ", H_total)
        println("Φ_total: ", Φ_total)

        if mode == :PROG

        elseif mode == :EKMAN

            if length(height_level_counts) != 3
                throw(ErrorException("Error: in mode :EKMAN, only three coarse layers are allowed."))
            end

        else
            throw(ErrorException("Error: Unknown mode :" * string(mode)))
        end

        s = zeros(Complex{Float64}, Nx, Ny)
        @. s = ϵ + gi.c_f * im


        return new(
            gi,
            Δt,
            Kh_barotropic,
            Kh_baroclinic,
            Nx, Ny, Nz, 
            z_bnd,
            H,
            H_total,
            Φ_total,
            mask,
            mode,
            ϵ,
            s,
        ) 
    end
    

end


