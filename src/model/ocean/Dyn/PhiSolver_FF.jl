using LinearAlgebra


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

        # Create coversion matrix and its inverse
        T_num = reshape(collect(1:length(mask2)), size(mask2)...)
        active_num_eff     = T_num[ mask2 .==1 ]

        eT_send_T = M.op.T_I_T[active_num_eff, :]
        
        dropzeros!(eT_send_T)

        T_send_eT = sparse(eT_send_T')
       
        # identity 
        eT_I_eT = eT_send_T * T_send_eT

        # Laplacian without boundary condition assigned 0 (borderfilter_T)
        T_Lap_no_bc_T   = M.filter_T * (M.T_DIVx_U * M.U_∂x_T + M.T_DIVy_V * M.V_∂y_T) * M.filter_T

        # Laplacian with boundary condition assigned 0 on F grid

        #println("size of T_Lap_no_bc_T: ", size(T_Lap_no_bc_T))
        #println("size of M.T_interp_F: ", size(M.T_interp_F))
        #println("size of M.borderfilter_F: ", size(M.borderfilter_F))
        #println("size of M.F_interp_T: ", size(M.F_interp_T))

        T_Lap_T = T_Lap_no_bc_T * M.T_interp_F * M.borderfilter_F * M.F_interp_T

        # ceT_Lap_ceT will get incorrect Laplacian since some info is lost during compression
        #eT_Lap_eT = eT_send_T * T_Lap_T * T_send_eT
        eT_Lap_eT = eT_send_T * T_Lap_T * T_send_eT

        # Modified Laplacian
        eT_MoLap_eT = eT_Lap_eT - α * eT_I_eT

        MoLap = lu(eT_MoLap_eT)
        tool_mtx = (
            MoLap = lu(eT_MoLap_eT),
        )
        
        return new(
            size(eT_MoLap_eT)[1],
            M,
            α,
            eT_send_T,
            T_send_eT,
            eT_Lap_eT,
            eT_MoLap_eT,
            T_Lap_T,
            tool_mtx,
        )

    end
end


