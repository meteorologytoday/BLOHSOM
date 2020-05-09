using LinearAlgebra


mutable struct PhiSolver

    eF_length    :: Int64

    M :: DynamicAdvSpeedUpMatrix
    α :: Float64
            
    eF_interp_T  :: AbstractArray{Float64, 2}
    eF_send_F    :: AbstractArray{Float64, 2}
    F_send_eF    :: AbstractArray{Float64, 2}
    eF_Lap_eF    :: AbstractArray{Float64, 2}
    eF_MoLap_eF  :: AbstractArray{Float64, 2}

    F_Lap_F
   
    tool_mtx

    wksp

    function PhiSolver(;
        gi             :: PolelikeCoordinate.GridInfo,
        mask2          :: AbstractArray{Float64, 2},
        α              :: Float64,  # ((Δt)^2 * H)^(-1)
        M              :: Union{DynamicAdvSpeedUpMatrix, Nothing} = nothing
    )
        if M == nothing
            M = DynamicAdvSpeedUpMatrix(;
                gi = gi,
                Nz = 1,
                mask2 = mask2,
            )
        end

        Nx = gi.Nx
        Ny = gi.Ny

         
        # Create coversion matrix and its inverse
        mask2_F = reshape( M.borderfilter_F * ones(Float64, Nx*(Ny+1)), Nx, Ny+1)
        F_num = reshape(collect(1:length(mask2_F)), size(mask2_F)...)
        num_eff = F_num[ mask2_F .==1 ]

        eF_send_F = M.op.F_I_F[num_eff, :]
       

        println("!!!!!!!!!!!!size: ", size(eF_send_F))
        println("!!!!!!!!!!!! SUM of mask2_F: ", sum(mask2_F))
 
        dropzeros!(eF_send_F)

        F_send_eF = sparse(eF_send_F')
       
        # identity 
        eF_I_eF = eF_send_F * F_send_eF


        eF_interp_T = eF_send_F * M.F_interp_T

        F_Lap_F = M.F_Lap_F * M.borderfilter_F

        # ceF_Lap_ceF will get incorrect Laplacian since some info is lost during compression
        #eF_Lap_eF = eF_send_T * T_Lap_T * T_send_eF
        eF_Lap_eF = eF_send_F * F_Lap_F * F_send_eF

        # Modified Laplacian
        eF_MoLap_eF = eF_Lap_eF - α * eF_I_eF


        tool_mtx = (
            MoLap = lu(eF_MoLap_eF),
            Lap = lu(eF_Lap_eF),
        )
        

        eF_length = size(eF_MoLap_eF)[1]

        wksp = (
            rhs_T = zeros(Float64, Nx, Ny),
            rhs_eF = zeros(Float64, eF_length),
            lhs_eF = zeros(Float64, eF_length),
            lhs_F  = zeros(Float64, Nx * (Ny+1)),
        )

        return new(
            eF_length,
            M,
            α,
            eF_interp_T,
            eF_send_F,
            F_send_eF,
            eF_Lap_eF,
            eF_MoLap_eF,
            F_Lap_F,
            tool_mtx,
            wksp,
        )

    end
end

function solvePhi!(
    solver :: PhiSolver,
    input  :: AbstractArray{Float64},  # Φ_aux
    output :: AbstractArray{Float64}, 
)

    wksp = solver.wksp

    @. wksp.rhs_T = - input * solver.α
    
    mul!( wksp.rhs_eF, solver.eF_interp_T, view(wksp.rhs_T, :))
    
    ldiv!(wksp.lhs_eF, solver.tool_mtx.Lap, wksp.rhs_eF)

    mul!(wksp.lhs_F,      solver.F_send_eF,    wksp.lhs_eF)
    mul!(view(output, :), solver.M.T_interp_F, wksp.lhs_F )

end
