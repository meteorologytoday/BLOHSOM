mutable struct DynModel

    env     :: DynEnv
    state   :: DynState
    forcing :: DynForcing
    core    :: DynCore
    
    
    function DynModel(;
        gi :: PolelikeCoordinate.GridInfo,
        Δt,
        Kh_barotropic,
        Kh_baroclinic,
        z_bnd,
        mode,
        mask=nothing,
    )

        env = DynEnv(;
            gi = gi,
            Δt = Δt,
            Kh_barotropic = Kh_barotropic,
            Kh_baroclinic = Kh_baroclinic,
            Nx = gi.Nx,
            Ny = gi.Ny,
            z_bnd = z_bnd,
            mode  = mode,
            mask=mask,
        )

        state = DynState(env)
        forcing = DynForcing(env)
        core  = DynCore(env, state)

        #=
        state.u_c .= 1000
        state.v_c .= 1000

        uu = copy(state.u_c)
        vv = copy(state.v_c)

        mul3!(uu, core.s_ops.filter_U, state.u_c)
        mul3!(vv, core.s_ops.filter_V, state.v_c)
        
        state.u_c .-= uu
        state.v_c .-= vv
       
        state.u_c[uu .== 0] .= NaN 
        state.v_c[vv .== 0] .= NaN 
        =#


        return new(
            env,
            state,
            forcing,
            core,
        )

    end
end


