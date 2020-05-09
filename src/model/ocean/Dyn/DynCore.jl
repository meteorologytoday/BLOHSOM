mutable struct Layers
    u :: Array{SubArray}
    v :: Array{SubArray}
    u_total :: Array{SubArray}
    v_total :: Array{SubArray}
    function Layers()
        return new()
    end
end



mutable struct DynCore    # Adam Bashford

    c_ops    :: Union{DynamicAdvSpeedUpMatrix, Nothing}
    s_ops    :: Union{DynamicAdvSpeedUpMatrix, Nothing}
    va       :: Union{VerticalAverager, Nothing}
    Φ_solver :: PhiSolver
    diffusion_solver_Kh_barotropic :: DiffusionSolver
    diffusion_solver_Kh_baroclinic :: DiffusionSolver

    G_idx    :: Dict
    wksp     :: Workspace

    u_aux    :: AbstractArray{Float64, 3}
    v_aux    :: AbstractArray{Float64, 3}
    Φ_aux    :: AbstractArray{Float64, 2}

    ∂B∂x     :: AbstractArray{Float64, 3}
    ∂B∂y     :: AbstractArray{Float64, 3}

    layers   :: Layers

    # Used by EKMAN scheme
    M̃    :: AbstractArray{Complex{Float64}, 2}
    ṽ_ek :: AbstractArray{Complex{Float64}, 2}
    ṽ_rf :: AbstractArray{Complex{Float64}, 2}

 
    function DynCore(env, state)
        
        Nx = env.Nx
        Ny = env.Ny
        Nz = env.Nz

        va = VerticalAverager(
            z_bnd_f = env.z_bnd,
            height_level_counts = ones(Int64, Nz) 
        )
        
        println("Making Spatial Operators")
        @time s_ops = DynamicAdvSpeedUpMatrix(;
                gi = env.gi,
                Nz = 1,
                mask2 = env.mask,
        )

        c_ops = s_ops

        Φ_solver = PhiSolver(
            gi    = env.gi,
            mask2 = env.mask,
            α     = 1.0 / (env.Δt^2 * env.Φ_total),
            M     = s_ops,
        )

        diffusion_solver_Kh_barotropic = DiffusionSolver(;
            gi = env.gi,
            M  = s_ops,
            K  = env.Kh_barotropic,
            Δt = env.Δt,
            mask2 = env.mask,
        )

        diffusion_solver_Kh_baroclinic = DiffusionSolver(;
            gi = env.gi,
            M  = s_ops,
            K  = env.Kh_baroclinic,
            Δt = env.Δt,
            mask2 = env.mask,
        )




        #         now   Δt-ago  2Δt-ago
        G_idx = Dict(
            :now           =>  1,
            :one_Δt_ago    =>  2,
            :two_Δt_ago    =>  3,
        )

        println("Creating Workspaces")
        wksp = Workspace(;
            Nx = Nx,
            Ny = Ny,
            Nz = Nz,
            T = 5,
            U = 6,
            V = 6,
            sT = 7,
            sU = 5,
            sV = 5,
            shape = :xyz,
        ) 

        u_aux = zeros(Float64, Nx, Ny,   Nz)
        v_aux = zeros(Float64, Nx, Ny+1, Nz)
        Φ_aux = zeros(Float64, Nx, Ny)

        ∂B∂x  = zeros(Float64, Nx, Ny,   Nz)
        ∂B∂y  = zeros(Float64, Nx, Ny+1, Nz)
        
        # making layer-wise views
        layers = Layers()
        for (var, ref) in Dict(
            :u => state.u,
            :v => state.v,
            :u_total => state.u_total,
            :v_total => state.v_total,
        )
            setfield!(layers, var, genLayerView(ref))
        end

        M̃    = zeros(Complex{Float64}, Nx, Ny)
        ṽ_ek = zeros(Complex{Float64}, Nx, Ny)
        ṽ_rf = zeros(Complex{Float64}, Nx, Ny)

        new(
            c_ops,
            s_ops,
            va,
            Φ_solver,
            diffusion_solver_Kh_barotropic,
            diffusion_solver_Kh_baroclinic,
            G_idx,
            wksp,
            u_aux,
            v_aux,
            Φ_aux,
            ∂B∂x,
            ∂B∂y,
            layers,
        )
    end
end

function genLayerView(
    arr :: AbstractArray{T, 3}
) where T

    # vectorization the for loop is not faster

    s = size(arr)

    _, _, Nz = size(arr)
    view_arr  = Array{SubArray}(undef, Nz)

    for k=1:Nz
        view_arr[k] = view(arr, :, :, k)
    end

    return view_arr
end


