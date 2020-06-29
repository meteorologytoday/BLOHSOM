
include("../../../share/PolelikeCoordinate.jl")

module Dyn

    using Formatting
    using LinearAlgebra    
    using ..PolelikeCoordinate
    using Statistics: mean

    include("../../../share/constants.jl")
    include("../../../share/ocean_state_function.jl")


    macro unpack(model)
        return esc(:( 
            co = $(model).core;
            st = $(model).state;
            fr = $(model).forcing;
            ev = $(model).env;
        ))
    end
 

    @inline function mul2!(
        a :: AbstractArray{Float64, 2},
        b :: AbstractArray{Float64, 2},
        c :: AbstractArray{Float64, 2},
    )

        mul!(view(a, :), b, view(c, :))

    end
 
    @inline function mul3!(
        a :: AbstractArray{Float64, 3},
        b :: AbstractArray{Float64, 2},
        c :: AbstractArray{Float64, 3},
    )

        for k=1:size(a)[3]
            
            mul!(
                view(view(a, :, :, k), :),
                b,
                view(view(c, :, :, k), :),
            )
        end

    end
 
    const α_AB = 1.0 /  2.0 
    const β_AB = 5.0 / 12.0  
    @inline function ABIII(a, aa, aaa)
        return ( 1 + α_AB + β_AB ) * a - (α_AB + 2*β_AB) * aa + β_AB * aaa
    end
   

    include("../MatrixOperators.jl")
    include("../VerticalAverager.jl")
    include("../Workspace.jl")
    include("../allocat.jl")

    include("AdvectionSpeedUpMatrix_dyn.jl")
    #include("PhiSolver.jl")
    #include("PhiSolver_F.jl")
    #include("PhiSolver.jl.old")
    include("PhiSolver_T_bc0.jl")
    include("DiffusionSolver.jl")



    include("DynEnv.jl")
    include("DynState.jl")
    include("DynForcing.jl")
    include("DynCore.jl")
    include("DynModel.jl")
    include("step_model.jl")


    include("varlist.jl")

   
    function stepModel!(
        m :: DynModel,
    )
        @unpack m

        reset!(co.wksp)

        if ev.mode == :PROG
            
            advectDynamic!(m)

        elseif env.mode == :EKMAN

            assignEkmanFlow!(m) 

        end

    end


end
