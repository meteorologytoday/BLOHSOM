mutable struct TmdModel

    env     :: TmdEnv
    state   :: TmdState
    core    :: TmdCore
    bridge_state :: Union{BridgeState, Nothing}

    function TmdModel(;
        env,
        state,
    )
        core  = TmdCore(env, state)

        if env.z_bnd_f != nothing
            bridge_state = BridgeState(env.z_bnd_f, env.height_level_counts)
        else
            bridge_state = nothing
        end

        return new(
            env, state, core, bridge_state,
        )

    end
end

