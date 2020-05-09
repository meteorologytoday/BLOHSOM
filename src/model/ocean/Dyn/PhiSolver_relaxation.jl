
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
        #mask2_exclude_bnd = reshape(M.borderfilter_T * view(mask2, :), Nx, Ny)

        # Create coversion matrix and its inverse
        T_num = reshape(collect(1:length(mask2)), size(mask2)...)
        active_num             = T_num[ mask2 .==1 ]
        #active_num_exclude_bnd = T_num[ mask2_exclude_bnd .==1 ]


        eT_send_T  = M.op.T_I_T[active_num, :]
        #ceT_send_T = M.op.T_I_T[active_num_exclude_bnd, :]
        
        dropzeros!(eT_send_T)
        #dropzeros!(ceT_send_T)

        T_send_eT  = sparse(eT_send_T')
        #T_send_ceT = sparse(ceT_send_T')
       
        # identity 
        #ceT_I_ceT = ceT_send_T * T_send_ceT
        eT_I_eT = eT_send_T * T_send_eT

        # Laplacian on T grids with gradient equals zero at boundaries (U, V grid boundaries)
        T_Lap_T   = M.T_DIVx_U * M.U_∂x_T + M.T_DIVy_V * M.V_∂y_T


        # ceT_Lap_ceT will get incorrect Laplacian since some info is lost during compression
        eT_Lap_eT = eT_send_T * T_Lap_T * T_send_eT

        dropzeros!(eT_Lap_eT)

        #println(Array(eT_Lap_eT))

        # Modified Laplacian
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

        MoLap = lu(eT_MoLap_eT)
        tool_mtx = (
        #    Lap   = lu(eT_Lap_eT),
            MoLap = lu(eT_MoLap_eT),
        )
        
        return new(
            size(eT_MoLap_eT)[1],
            M,
            α,
        #    ceT_send_T,
            eT_send_T,
        #    T_send_ceT,
            T_send_eT,
            eT_Lap_eT,
            eT_MoLap_eT,
            T_Lap_T,
            tool_mtx,
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
    mul!( wksp.rhs_eT, solver.eT_send_T, view(wksp.rhs_T, :))

    for k=1:inf
        mul!(wksp.lhs_F, solver.F_send_eF,    wksp.lhs_eF)
        
        
    end

    
    ldiv!(wksp.lhs_eF, solver.tool_mtx.MoLap, wksp.rhs_eF)

    mul!(wksp.lhs_F,      solver.F_send_eF,    wksp.lhs_eF)
    mul!(view(output, :), solver.M.T_interp_F, wksp.lhs_F )

end