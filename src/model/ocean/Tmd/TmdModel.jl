mutable struct TmdModel

    env     :: TmdEnv
    state   :: TmdState
    core    :: TmdCore

    function TmdModel(;
        env,
        state,
    )
        core  = TmdCore(env, state)
        
        return new(
            env, state, core,
        )

    end
end

