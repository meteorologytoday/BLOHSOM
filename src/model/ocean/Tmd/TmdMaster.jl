mutable struct TmdMaster

    pplan   :: ParallelPlan
    env     :: TmdEnv
    state   :: TmdState
    current_substep :: Int64

    function TmdMaster(;
        gf,
        pplan,
        Δt,
        substeps,
        z_bnd,
        topo        = nothing,
        Kh_X        = [1e3,  1e3 ],
        Kv_X        = [1e-5, 1e-5],
        we_max      = 1e-2,
        R           = 0.58,  # See Paulson and Simpson (1977) Type I clear water
        ζ           = 23.0,  # See Paulson and Simpson (1977) Type I clear water
        MLT_rng     = [10.0, 1000.0],
        NX_passive  = 0,
        t_X_wr      = nothing,
        X_wr        = nothing,
        mask2       = nothing,
        MLT_scheme  = :dynamic,
        radiation_scheme = :exponential_decay,
        convective_adjustment = true,
        use_Qflux      = true,
        finding_Qflux  = false,
        z_bnd_f = nothing,
        height_level_counts = nothing,
    )

        # create master state
        gi = PolelikeCoordinate.genGridInfo(gf)

        master_env = TmdEnv(;
            gi         = gi,
            update_yrng_T   = nothing,
            Δt         = Δt,
            substeps   = substeps,
            z_bnd      = z_bnd,
            topo       = topo,
            Kh_X       = Kh_X,
            Kv_X       = Kv_X,
            we_max     = we_max,
            R          = R,
            ζ          = ζ,
            MLT_rng    = MLT_rng, 
            NX_passive = NX_passive,
            t_X_wr     = t_X_wr,
            X_wr       = X_wr,
            mask2      = mask2,
            MLT_scheme = MLT_scheme,
            radiation_scheme = radiation_scheme,
            convective_adjustment = convective_adjustment,
            use_Qflux  = use_Qflux,
            finding_Qflux = finding_Qflux,
            z_bnd_f = z_bnd_f,
            height_level_counts = height_level_counts,
        )

        master_state = TmdState(master_env)

        @sync for (i, p) in enumerate(pplan.pids)
            
            visible_yrng_T = pplan.visible_yrngs_T[i]
            visible_yrng_V = pplan.visible_yrngs_V[i]
            update_yrng_T  = pplan.update_yrngs_T[i]

            slave_gi = PolelikeCoordinate.genGridInfo(
                gf,
                sub_yrng=visible_yrng_T,
            )
 
            slave_env = TmdEnv(;
                gi         = slave_gi,
                update_yrng_T = update_yrng_T,
                Δt         = Δt,
                substeps   = substeps,
                z_bnd      = z_bnd,
                topo       = ( topo != nothing) ? topo[:, visible_yrng_T] : topo,
                Kh_X       = Kh_X,
                Kv_X       = Kv_X,
                we_max     = we_max,
                R          = R,
                ζ          = ζ,
                MLT_rng    = MLT_rng, 
                NX_passive = NX_passive,
                t_X_wr     = t_X_wr,
                X_wr       = ( X_wr != nothing ) ? X_wr[:, :, :, visible_yrng_T] : X_wr,
                mask2      = mask2[:, visible_yrng_T],
                MLT_scheme = MLT_scheme,
                radiation_scheme = radiation_scheme,
                convective_adjustment = convective_adjustment,
                use_Qflux  = use_Qflux,
                finding_Qflux = finding_Qflux,
                z_bnd_f = z_bnd_f,
                height_level_counts = height_level_counts,
            )
   
            slave_state = createMirror(
                master_state,
                master_env.Ny,
                visible_yrng_T,
                visible_yrng_V,
            ) 


            @spawnat p let
                createSlaveTmdModel(;
                    env   = slave_env,
                    state = slave_state,
                )

            end 
 
        end
        
        return new(
            pplan, master_env, master_state, 1,
        )

    end
end

function createSlaveTmdModel(;
    env :: TmdEnv,
    state :: TmdState,
)
                
    println("### Define tmd_model at pid: ", myid())
    global tmd_model = TmdModel(;
        env = env,
        state = state,
    )

end


function setupBridgeState!(
    B_to_dyn :: AbstractArray{Float64, 3},
    uc_from_dyn :: AbstractArray{Float64, 3},
    vc_from_dyn :: AbstractArray{Float64, 3},
)
    global tmd_model
    update_yrng_T = tmd_model.env.update_yrng_T

    # received variables are in :xyz
    tmd_model.bridge_state.uc_from_dyn = view(uc_from_dyn, :, update_yrng_T, :)
    tmd_model.bridge_state.vc_from_dyn = view(vc_from_dyn, :, update_yrng_T, :)
    tmd_model.bridge_state.B_to_dyn = view(B_to_dyn, :, update_yrng_T, :)

end

