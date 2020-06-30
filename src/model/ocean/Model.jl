mutable struct Model

    env            :: OcnEnv
    tmd_pplan      :: Tmd.ParallelPlan
    dyn_engine     :: Dyn.DynModel
    tmd_engine     :: Tmd.TmdMaster
#    data_table     :: DataTable

    function Model(
        env :: OcnEnv,
    )
         
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
        tmd_pplan      = Tmd.ParallelPlan(env.Ny)
        tmd_engine = Tmd.TmdMaster(
            gf       = env.hrgrid,
            pplan    = tmd_pplan, 
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
            NX_passive = env.NX_passive,
            t_X_wr     = env.t_X_wr,
            X_wr       = env.X_wr,
            MLT_scheme = env.MLT_scheme,
            radiation_scheme = Symbol(env.radiation_scheme),
            convective_adjustment = env.convective_adjustment,
            use_Qflux     = env.use_Qflux,
            finding_Qflux = env.finding_Qflux,
            z_bnd_f = env.z_bnd_f,
            height_level_counts = env.height_level_counts,
        )         

        @sync for p in tmd_pplan.pids
                @spawnat p let 
                    Tmd.setupBridgeState!(
                        dyn_engine.state.B, 
                        dyn_engine.state.u_total,
                        dyn_engine.state.v_total,
                    )
                end
        end

        #=
        # Registering variables
        descs_X = (
            (:X,     :fT, :zxy, tmd_engine.state.X,    Float64),
            (:X_ML,  :sT, :xy,  tmd_engine.state.X_ML, Float64),
        )
     
        descs_noX = (

            (:h_ML,  :sT, :xy,  tmd_engine.state.h_ML, Float64),
            (:b_ML,  :sT, :xy,  tmd_engine.state.b_ML, Float64),
            (:b   ,  :fT, :zxy, tmd_engine.state.b,    Float64),
            (:B   ,  :fT, :zxy, tmd_engine.state.B     Float64),
            (:FLDO,  :sT, :xy,  tmd_engine.state.FLDO,   Int64),
            
            # These are used by dyn_core 
            (:u_c, :cU, :xyz, tmdFloat64),
            (:v_c, :cV, :xyz, Float64),
            (:B_c,       :cT, :xyz, Float64),
            (:Φ,         :sT, :xy,  Float64),
            (:∂B∂x,      :cU, :xyz, Float64),
            (:∂B∂y,      :cV, :xyz, Float64),

            # Forcings and return fluxes to coupler
            (:SWFLX,   :sT, :xy,  Float64),
            (:NSWFLX,  :sT, :xy,  Float64),
            (:TAUX,    :sT, :xy,  Float64),
            (:TAUY,    :sT, :xy,  Float64),
            (:IFRAC,   :sT, :xy,  Float64),
            (:VSFLX,   :sT, :xy,  Float64),
            (:QFLX_T,  :sT, :xy,  Float64),
            (:QFLX_S,  :sT, :xy,  Float64),
            (:T_CLIM,  :sT, :xy,  Float64),
            (:S_CLIM,  :sT, :xy,  Float64),
            (:MLT,     :sT, :xy,  Float64),
            (:w_W,     :fW, :xyz,  Float64),
        ) 


        for (id, grid, shape, dtype) in descs_X
            regVariable!(model.shared_data, model.env, id, grid, shape, dtype, has_Xdim=true)
        end

        for (id, grid, shape, dtype) in descs_noX
            regVariable!(model.shared_data, model.env, id, grid, shape, dtype, has_Xdim=false)
        end



        =#

        return new(
            env,
            tmd_pplan,
            dyn_engine,
            tmd_engine,
        )
    end
end

