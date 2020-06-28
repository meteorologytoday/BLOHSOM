mutable struct TmdModel

    pplan   :: ParallelPlan
    env     :: TmdEnv
    core    :: TmdCore
    state   :: Array{TmdState, 1}
    diag    :: TmdDiag
    forcing :: TmdForcing

    function TmdModel(;
        workers,
        gi,
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


        env = TmdEnv(;
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

        pplan = ParallelPlan(env, workers)
        states  = [TmdState(env) for _ in 1:2]
        diag    = TmdDiag(env)
        forcing = TmdForcing(env)
 
        # Create slaves based on plan
        slaves = [ TmdSlave(
            pplan,
            states,
            diag,
            forcing
        ) for _ in 1:length(pplan.pids) ]



      
        # core now should belong to each slave 
        # core  = TmdCore(env, state, diag, forcing)



        return new(
            env, state, diag, forcing
        )

    end
end
