mutable struct Model

    env            :: OcnEnv
    dyn_engine     :: DynModel
    tmd_engine     :: TmdMaster

    job_dist_info  :: JobDistributionInfo

    function Model(
        env :: OcnEnv,
    )
        
        tmd_pplan      = Tmd.ParallelPlan(env.Ny)

        # setup dyn engine
        gi = PolelikeCoordinate.genGridInfo(
            env.hrgrid,
        )

        dyn_engine = Dyn.DynModel(
            mode                    = env.flow_scheme,
            gi                      = gi,
            Δt                      = env.Δt / env.substeps_dyn,
            Kh_barotropic           = env.Kh_m_barotropic,
            Kh_baroclinic           = env.Kh_m_baroclinic,
            Kv_baroclinic           = env.Kv_m_baroclinic,
            τ_barotropic_bottomfric = env.τ_barotropic_bottomfric,
            τ_barotropic_coastfric  = env.τ_barotropic_coastfric,
            z_bnd                   = env.z_bnd_c,
            mask                    = env.mask2_deep,
        )         



        # setup tmd engine
        tmd_engine = Tmd.TmdMaster(
            gi       = gi,
            Δt       = env.Δt,
            substeps = env.substeps_tmd,
            z_bnd    = env.z_bnd_f,
            topo     = env.topo,
            mask2    = env.mask2,
            Kh_X     = env.Kh_X,
            Kv_X     = env.Kv_X,
            we_max   = env.we_max,
            R        = env.R,
            ζ        = env.ζ,
            MLT_rng  = env.MLT_rng,
            NX_passive = env. NX_passive,
            t_X_wr     = env.t_X_wr,
            X_wr       = X_wr,
            MLT_scheme = env.MLT_scheme,
            radiation_scheme = Symbol(env.radiation_scheme),
            convective_adjustment = env.convective_adjustment,
            use_Qflux     = env.use_Qflux,
            finding_Qflux = env.finding_Qflux,
        )         

        # Register data exchange plan
        data_exchanger = DataExchanger([
           :TMD2DYN, :DYN2TMD, 
        ])

        hasXdim = true
        noXdim  = false

        bindings = (
            ((:u_c,       :cU, :xyz, bd[:u_c], noXdim), :u_c),
            ((:v_c,       :cV, :xyz, bd[:v_c], noXdim), :v_c),
            ((:B_c,       :cT, :xyz, bd[:B_c],       noXdim), :B_c),
            
            # T, S, FLDO and such
            ((:X,    :fT, :zxy, s.X,    hasXdim), :X   ),
            ((:X_ML, :sT, :xy , s.X_ML, hasXdim), :X_ML),
            ((:h_ML, :sT, :xy , s.h_ML, noXdim),  :h_ML),
            ((:FLDO, :sT, :xy , s.FLDO, noXdim),  :FLDO),
            ((:b,    :fT, :zxy, s.b,    noXdim),  :b   ),
            ((:b_ML, :sT, :xy,  s.b_ML, noXdim),  :b_ML),
            ((:B,    :fT, :zxy, s.B,    noXdim),  :B   ),

            # forcings
            ((:SWFLX,  :sT, :xy,  f.swflx,  noXdim), :SWFLX),
            ((:NSWFLX, :sT, :xy,  f.nswflx, noXdim), :NSWFLX),
            ((:w_W,    :fW, :zxy, f.w_W, noXdim),    :w_W),

        )

        for (here_args, there_key) in bindings
           
            here = DataUnit(here_args...)
            println("Doing : ", here.id, "; ", du_there[there_key].id, "; size: ", size(here.data)) 
            
            here_pull_yrng  = Colon()

            if here.grid in (:fV, :cV, :sV)
                there_pull_yrng = ysi.pull_fr_rng[1]:(ysi.pull_fr_rng[end]+1)
                here_push_yrng  = ysi.push_fr_rng[1]:(ysi.push_fr_rng[end]+1)
                there_push_yrng = ysi.push_to_rng[1]:(ysi.push_to_rng[end]+1)
            else
                there_pull_yrng  = ysi.pull_fr_rng
                here_push_yrng   = ysi.push_fr_rng
                there_push_yrng  = ysi.push_to_rng
            end

            createBinding!(
                de,
                there_key,
                here,
                du_there[there_key],
                here_pull_yrng,
                there_pull_yrng,
                here_push_yrng,
                there_push_yrng,
            )
        end

        #=
        for (label, du_ids) in groups
            for id in du_ids
                addToGroup!(de, id, label)
            end
        end
        =#

        #println("done.")
    end



    end
end

function regVariable!(
    sd       :: SharedData,
    env      :: OcnEnv,
    id       :: Symbol,
    grid     :: Symbol,
    shape    :: Symbol,
    dtype    :: DataType;
    data     :: Union{Nothing, SharedArray{T}} = nothing,
    has_Xdim :: Bool = false,
) where T


    env = sd.env
    Nx, Ny, Nz_f, Nz_c = env.Nx, env.Ny, env.Nz_f, env.Nz_c
    NX = env.NX


    if haskey(sd.data_units, id)
        throw(ErrorException("Error: variable id " * String(id) *  " already exists."))
    end

    dim = Dict(
        :fT => [Nx  , Ny  , Nz_f  ],
        :fU => [Nx  , Ny  , Nz_f  ],
        :fV => [Nx  , Ny+1, Nz_f  ],
        :fW => [Nx  , Ny  , Nz_f+1],
        :cT => [Nx  , Ny  , Nz_c  ],
        :cU => [Nx  , Ny  , Nz_c  ],
        :cV => [Nx  , Ny+1, Nz_c  ],
        :cW => [Nx  , Ny  , Nz_c+1],
        :sT => [Nx  , Ny  ],
        :sU => [Nx  , Ny  ],
        :sV => [Nx  , Ny+1],
    )[grid]

    if ! (shape in (:xy, :xyz, :zxy))
        throw(ErrorException("Error: only :xy, :xyz, :zxy are accepted in shape"))
    end

    if grid in (:sT, :sU, :sV)
        if shape in (:xyz, :zxy)
            throw(ErrorException("Error: grid and shape mismatch. Slab grid should only be with :xy shape."))
        end
    else
        if shape == :zxy
            dim = circshift(dim, 1)
        end
    end

    if has_Xdim
        push!(dim, NX)
    end

    if ! (dtype in (Float64, Int64))
        throw(ErrorException("Invalid data type. Only Float64 and Int64 are accepted"))
    end

    if data == nothing
        data = SharedArray{dtype}(dim...)
    else

        if Tuple(dim) != size(data)
            println("Expect ", dim)
            println("Get ", size(data))
            throw(ErrorException("Provided data does not have correct dimension."))
        end

        if dtype != T
            throw(ErrorException("dtype and provided data does not match."))
        end
    end

    sd.data_units[id] = DataUnit(
        id,
        grid,
        shape,
        data,
        has_Xdim,
    )

    sd.flags[id] = 0
end


function setupBinding!(
    slave :: TmdSlave,
)

    #println("TmdSlave setupBinding...")

    de = slave.data_exchanger
    sd = slave.shared_data
    m  = slave.model
    du_there = sd.data_units

    bd = slave.buffer_data
    s  = m.state
    f  = m.forcing

    ysi = slave.y_split_info

    hasXdim = true
    noXdim  = false

    bindings = (
        ((:u_c,       :cU, :xyz, bd[:u_c], noXdim), :u_c),
        ((:v_c,       :cV, :xyz, bd[:v_c], noXdim), :v_c),
        ((:B_c,       :cT, :xyz, bd[:B_c],       noXdim), :B_c),
        
        # T, S, FLDO and such
        ((:X,    :fT, :zxy, s.X,    hasXdim), :X   ),
        ((:X_ML, :sT, :xy , s.X_ML, hasXdim), :X_ML),
        ((:h_ML, :sT, :xy , s.h_ML, noXdim),  :h_ML),
        ((:FLDO, :sT, :xy , s.FLDO, noXdim),  :FLDO),
        ((:b,    :fT, :zxy, s.b,    noXdim),  :b   ),
        ((:b_ML, :sT, :xy,  s.b_ML, noXdim),  :b_ML),
        ((:B,    :fT, :zxy, s.B,    noXdim),  :B   ),

        # forcings
        ((:SWFLX,  :sT, :xy,  f.swflx,  noXdim), :SWFLX),
        ((:NSWFLX, :sT, :xy,  f.nswflx, noXdim), :NSWFLX),
        ((:w_W,    :fW, :zxy, f.w_W, noXdim),    :w_W),

    )

    #=
    groups = Dict(
        :FR_DYN => (:u_c, :v_c,),
        :TO_DYN => (:b_c,),
        :BND    => (:X, :X_ML,),
        :FR_MAS => (:SWFLX, :NSWFLX),
        :TO_MAS => (:X, :X_ML, :h_ML, :FLDO),
    )
    =#
    #println("createBinding..")

    for (here_args, there_key) in bindings
       
        here = DataUnit(here_args...)
        println("Doing : ", here.id, "; ", du_there[there_key].id, "; size: ", size(here.data)) 
        
        here_pull_yrng  = Colon()

        if here.grid in (:fV, :cV, :sV)
            there_pull_yrng = ysi.pull_fr_rng[1]:(ysi.pull_fr_rng[end]+1)
            here_push_yrng  = ysi.push_fr_rng[1]:(ysi.push_fr_rng[end]+1)
            there_push_yrng = ysi.push_to_rng[1]:(ysi.push_to_rng[end]+1)
        else
            there_pull_yrng  = ysi.pull_fr_rng
            here_push_yrng   = ysi.push_fr_rng
            there_push_yrng  = ysi.push_to_rng
        end

        createBinding!(
            de,
            there_key,
            here,
            du_there[there_key],
            here_pull_yrng,
            there_pull_yrng,
            here_push_yrng,
            there_push_yrng,
        )
    end

    #=
    for (label, du_ids) in groups
        for id in du_ids
            addToGroup!(de, id, label)
        end
    end
    =#

    #println("done.")
end


