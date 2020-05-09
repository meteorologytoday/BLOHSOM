
mutable struct PhiSolver

    eT_length    :: Int64

    M :: DynamicAdvSpeedUpMatrix
    α :: Float64

    eT_send_T    :: AbstractArray{Float64, 2}
    T_send_eT    :: AbstractArray{Float64, 2}
    eT_Lap_eT    :: AbstractArray{Float64, 2}
    eT_MoLap_eT  :: AbstractArray{Float64, 2}
   
    T_Lap_T
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

        # need a mask excluding bnd points
        mask2_eff = reshape(M.borderfilter_T * view(mask2, :), Nx, Ny)

        # Create coversion matrix and its inverse
        T_num = reshape(collect(1:length(mask2_eff)), size(mask2_eff)...)
        active_num_eff     = T_num[ mask2_eff .==1 ]


        #eT_send_T  = M.op.T_I_T[active_num, :]
        eT_send_T = M.op.T_I_T[active_num_eff, :]
        
        #dropzeros!(eT_send_T)
        dropzeros!(eT_send_T)

        #T_send_eT  = sparse(eT_send_T')
        T_send_eT = sparse(eT_send_T')
       
        # identity 
        #ceT_I_ceT = ceT_send_T * T_send_ceT
        #eT_I_eT = eT_send_T * T_send_eT
        eT_I_eT = eT_send_T * T_send_eT

        # Laplacian on T grids with gradient equals zero at boundaries (U, V grid boundaries)
        #T_Lap_T   = M.T_DIVx_U * M.U_∂x_T + M.T_DIVy_V * M.V_∂y_T

        # Laplacian with boundary condition assigned 0 (borderfilter_T)
        T_Lap_T   = (M.T_DIVx_U * M.U_∂x_T + M.T_DIVy_V * M.V_∂y_T) * M.borderfilter_T


        # ceT_Lap_ceT will get incorrect Laplacian since some info is lost during compression
        #eT_Lap_eT = eT_send_T * T_Lap_T * T_send_eT
        eT_Lap_eT = eT_send_T * T_Lap_T * T_send_eT

#        dropzeros!(eT_Lap_eT)

        #println(Array(eT_Lap_eT))

        # Modified Laplacian
        #eT_MoLap_eT = eT_Lap_eT - α * eT_I_eT

        eT_MoLap_eT = eT_Lap_eT - α * eT_I_eT


#=
        # Testing       
        rr, cc = size(eT_Lap_eT)

#        println(eT_Lap_eT.nzval)

        r_cnt = 0
        for i=1:rr
            if all(eT_Lap_eT[i, :] .<= 1e-13)
                r_cnt += 1
            end
        end

        c_cnt = 0
        for j=1:cc
            if all(eT_Lap_eT[:, j] .<= 1e-13)
                c_cnt += 1
            end
        end

        println("Empty rows: ", r_cnt)
        println("Empty cols: ", c_cnt)
=#

        #MoLap = lu(eT_MoLap_eT)
        tool_mtx = (
        #    Lap   = lu(eT_Lap_eT),
            MoLap = lu(eT_MoLap_eT),
        )
 
        eT_length = size(eT_MoLap_eT)[1]
        wksp = (
            rhs_T = zeros(Float64, Nx, Ny),
            rhs_eT = zeros(Float64, eT_length),
            lhs_eT = zeros(Float64, eT_length),
            lhs_T  = zeros(Float64, Nx * Ny),
        )
       
        return new(
            eT_length,
            M,
            α,
            eT_send_T,
            T_send_eT,
            eT_Lap_eT,
            eT_MoLap_eT,
            T_Lap_T,
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
    output .= 0.0    
    
    mul!(wksp.rhs_eT, solver.eT_send_T, view(wksp.rhs_T,:))
    ldiv!(wksp.lhs_eT, solver.tool_mtx.MoLap, wksp.rhs_eT)

    mul!(view(output, :),      solver.T_send_eT,    wksp.lhs_eT)

end
