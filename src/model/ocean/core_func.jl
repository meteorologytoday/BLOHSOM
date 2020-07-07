function init!(
    ocn_env = Union{String, OceanEnv};
    snapshot :: Union{String, Nothing} = nothing,
)

    if typeof(ocn_env) == String
        ocn_env = loadOcnEnv(ocn_env)
    end

    model = Model(ocn_env)


    return model 

end


function stepModel!(
    model :: Model,
    write_restart :: Bool,
)

    env = model.env

    Dyn.stepModel!(model.dyn_engine)

    ##### tmd slave should distribute u,v to fine grids here #####
    @sync for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
        @spawnat pid projVelocity!(tmd_slave)
    end


    # this involves passing tracer through boundaries
    # so need to sync every time after it evolves
    for t=1:env.substeps_tmd

        @sync for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
            @spawnat pid stepModel!(tmd_slave)
        end

        touchTmd!(model, :TMDBND, :S2M)
        touchTmd!(model, :TMDBND, :M2S)

    end
    
    ##### tmd slave should calcaulte b of coarse grid #####
    @sync for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
        @spawnat pid let
            BLOHSOM.calCoarseBuoyancyPressure!(tmd_slave)
        end
    end

    @sync let
        @async touchDyn!(model, :END_DYN2MAS, :S2M)
        @async touchTmd!(model, :END_TMD2MAS, :S2M)
    end

    touchTmd!(model, :TMD2DYN, :S2M)

    #= 
    if write_restart
        writeRestart(
            dyn_slave,
            tcr_slave,
            mld_slave, 
        )
    
    =#
end


#=

function init!(
    ocn_env = Union{String, OceanEnv};
    snapshot :: Union{String, Nothing} = nothing,
)

    # TODO
    ocn_env = (typeof(ocn_env) == String) ? loadOcnEnv(env) : ocn_env
    
    shared_data    = SharedData(ocn_env)
    model = Model(ocn_env)



    println("Register Shared Data") 
    regSharedData!(model)
    
    # Potential redundant: seems like shared_data is already doing this
    #ocn_state      = OcnState(shared_data)


    println("Creating slaves on nodes")
    @sync let

        @spawnat job_dist_info.dyn_slave_pid let
            global dyn_slave = BLOHSOM.DynSlave(ocn_env, shared_data)

            BLOHSOM.setupBinding!(dyn_slave)

        end

        for (p, pid) in enumerate(job_dist_info.tmd_slave_pids)
            @spawnat pid let
                global tmd_slave = BLOHSOM.TmdSlave(
                    ocn_env,
                    shared_data,
                    job_dist_info.y_split_infos[p],
                )

                BLOHSOM.setupBinding!(tmd_slave)

            end
        end
    end

    regBindingGroups!(model)

    if snapshot != nothing

        snapshot_dyn = format("{:s}.dyn.nc", snapshot)
        snapshot_tmd = format("{:s}.tmd.nc", snapshot)
        
        Dyn.loadSnapshot!(snapshot_dyn)
        Tmd.loadSnapshot!(snapshot_tmd)

    end



    println("Slave created and data exchanger is set.")
    println("TODO: Skip read restart file for now.")

    #=
    if restart_file != nothing
        loadRestart(restart_file, shared_data)
    end
    =#

    #=
    @sync let 
        touchDyn!(model, :END_DYN2MAS, :S2M)
        touchTmd!(model, :END_TMD2MAS, :S2M)
    end

    @sync touchTmd!(model, :TMD2DYN, :S2M)
    =#
    return model

end


function regBindingGroups!(model::Model)
    groups = Dict(
        :DYN2TMD     => (:u_c, :v_c),
        :TMD2DYN     => (:B_c,),
        :TMDBND      => (:X, :X_ML, :h_ML, :FLDO),
        :END_DYN2MAS => (:Φ, :∂B∂x, :∂B∂y),
        :END_TMD2MAS => (:X, :X_ML, :h_ML, :FLDO, :b, :b_ML, :B, :w_W),
        :FORCING_MAS2DYN => (),
        :FORCING_MAS2TMD => (:SWFLX, :NSWFLX),
        :TEST_MAS2TMD => (:X, :X_ML),
    ) 

    @sync let
        for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
            @spawnat pid let 

                for group_label in [ :DYN2TMD, :TMD2DYN, :TMDBND, :END_TMD2MAS, :FORCING_MAS2TMD]
                    for binding_id in groups[group_label]
                        addToGroup!(tmd_slave.data_exchanger, binding_id, group_label)
                    end
                end
            end
        end
        
        @spawnat model.job_dist_info.dyn_slave_pid let 
            for group_label in [:DYN2TMD, :TMD2DYN, :END_DYN2MAS, :FORCING_MAS2DYN]
                for binding_id in groups[group_label]
                    addToGroup!(dyn_slave.data_exchanger, binding_id, group_label)
                end
            end
        end

    end
end

function regSharedData!(model::Model)

    descs_X = (
        (:X,     :fT, :zxy, Float64),
        (:X_ML,  :sT, :xy,  Float64),
    )
 
    descs_noX = (

        (:h_ML,  :sT, :xy,  Float64),
        (:b_ML,  :sT, :xy,  Float64),
        (:b   ,  :fT, :zxy, Float64),
        (:B   ,  :fT, :zxy, Float64),
        (:FLDO,  :sT, :xy,    Int64),
        
        # These are used by dyn_core 
        (:u_c, :cU, :xyz, Float64),
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

end


push_pull_relation = Dict(
    :S2M => :PUSH,
    :M2S => :PULL,
)

function touchTmd!(
    model :: Model,
    group_label :: Symbol,
    direction :: Symbol,
)

    direction = push_pull_relation[direction]

    
    @sync for (p, pid) in enumerate(model.job_dist_info.tmd_slave_pids)
        @spawnat pid let
            BLOHSOM.syncData!(tmd_slave.data_exchanger, group_label, direction)
        end
    end
end

function touchDyn!(
    model:: Model,
    group_label :: Symbol,
    direction :: Symbol,
)

    direction = push_pull_relation[direction]

    @sync @spawnat model.job_dist_info.dyn_slave_pid let
        BLOHSOM.syncData!(dyn_slave.data_exchanger, group_label, direction)
    end
end

function loadData!()
end
=#
