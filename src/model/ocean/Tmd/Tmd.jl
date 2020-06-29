module Tmd

    @inline function cyc(i::Int64, N::Int64)
        return mod(i-1, N) + 1
    end

    @inline function mul_autoflat!(
        a :: AbstractArray{Float64},
        b :: AbstractArray{Float64, 2},
        c :: AbstractArray{Float64},
    )
        mul!(view(a, :), b, view(c, :))
    end
 
    macro fast_extract(model)
        return esc(:( 
            co = $(model).core;
            st = $(model).state;
            ev = $(model).env;
        ))
    end

    macro loop_hor(model, idx1, idx2, stmts)
        return :( for grid_idx in 1:size($(esc(model)).core.valid_idx)[2]

            $(esc(idx1)) = $(esc(model)).core.valid_idx[1, grid_idx]
            $(esc(idx2)) = $(esc(model)).core.valid_idx[2, grid_idx]
            $(esc(stmts))

        end )
    end

    using SharedArrays
    using Distributed
    using SparseArrays
    using Formatting
    using LinearAlgebra    
    using ..PolelikeCoordinate
    using ..GridFiles
    using Statistics: mean

    include("../../../share/constants.jl")
    include("../../../share/ocean_state_function.jl")

    include("ParallelPlan.jl")
    include("../MatrixOperators.jl")
    include("../Workspace.jl")
    include("allocate.jl")
    include("AdvectionSpeedUpMatrix.jl")
    include("AccumulativeVariables.jl")
    include("TmdEnv.jl")
    include("TmdState.jl")
    include("TmdCore.jl")
    include("TmdModel.jl")
    include("TmdMaster.jl")

    # functions
    include("latent_heat_release_of_freezing.jl")
    include("columnwise_budget.jl")
    include("trivial_functions.jl")

    include("mld_calculation.jl")
    include("convective_adjustment.jl")
    include("diffusion.jl")
    include("mixUnmix.jl")
    include("calFLDOPartition.jl")
    include("columnwise_integration.jl")
    #include("deep_ocn_correction.jl")
    include("shortwave_radiation.jl")
    #include("flx_correction.jl")
    include("initialization.jl")
    include("set_ocean_column.jl")
    include("buoyancy.jl")
 
    include("varlist.jl")

    include("leonard1979.jl")
    include("step_model.jl")
    #include("step_model.jl.old")
    include("step_tmd_mixed_layer.jl")

    function stepModel!(
        m :: TmdMaster,
    )

        @sync let
            for p in m.pplan.pids
                @spawnat p let
                    reset!(tmd_model.core.wksp)
                end
            end

            if m.current_substep == 1
                @sync for p in m.pplan.pids
                    @spawnat p let
                        determineVelocity!(tmd_model)
                    end
                end
            end
 
            for p in m.pplan.pids
                @spawnat p let
                    advectTracer_part1!(tmd_model)
                end
            end
               
                #println("sum of nswflx: ", sum(fr.nswflx))
                #println("co.current_substep: ", co.current_substep)
        end

        @sync let
            for p in m.pplan.pids
                @spawnat p let
                    advectTracer_part2!(tmd_model)
                    #doMixedLayerDynamics!(tmd_model)
                end
            end
        end


        if m.current_substep == m.env.substeps
            @sync let
                for p in m.pplan.pids
                    @spawnat p let
                        # do slow processes
                        #calBuoyancyPressure!(tmd_model)
                    end
                end
            end
        end

        if m.current_substep != m.env.substeps
            flag = :INTER_STEP
            m.current_substep += 1
        else
           flag = :FINAL_STEP
           m.current_substep = 1
        end
      
        return flag 
    end

end
