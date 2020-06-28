mutable struct TmdMaster

    env     :: TmdEnv
    state   :: TmdState

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
    )

        # create master state
        gi = PolelikeCoordinate.genGridInfo(gf)

        master_env = TmdEnv(;
            gi         = gi,
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
        )

        master_state = TmdState(master_env)

        for (i, p) in enumerate(pplan.pids)

            vis_rng_T = pplan.visible_rngs_T[i]
            vis_rng_V = pplan.visible_rngs_V[i]

            slave_gi = PolelikeCoordinate.genGridInfo(
                gf,
                sub_yrng=vis_rng_T,
            )
 
            slave_env = TmdEnv(;
                gi         = slave_gi,
                Δt         = Δt,
                substeps   = substeps,
                z_bnd      = z_bnd,
                topo       = ( topo != nothing) ? topo[:, vis_rng_T] : topo,
                Kh_X       = Kh_X,
                Kv_X       = Kv_X,
                we_max     = we_max,
                R          = R,
                ζ          = ζ,
                MLT_rng    = MLT_rng, 
                NX_passive = NX_passive,
                t_X_wr     = t_X_wr,
                X_wr       = ( X_wr != nothing ) ? X_wr[:, :, :, vis_rng_T] : X_wr,
                mask2      = mask2[:, vis_rng_T],
                MLT_scheme = MLT_scheme,
                radiation_scheme = radiation_scheme,
                convective_adjustment = convective_adjustment,
                use_Qflux  = use_Qflux,
                finding_Qflux = finding_Qflux,
            )
   
            slave_state = createMirror(
                master_state,
                master_env.Ny,
                vis_rng_T,
                vis_rng_V,
            ) 

            @spawnat p let
                global tmd_model = TmdModel(;
                    env   = slave_env,
                    state = slave_state,
                )
            end 
 
        end
        
        return new(
            master_env, master_state,
        )

    end
end
