
using LinearAlgebra

mutable struct DiffusionSolver

    eT_length    :: Int64

    M :: DynamicAdvSpeedUpMatrix
   
    K  :: Float64
    Δt :: Float64
    α  :: Float64

    eT_send_T    :: AbstractArray{Float64, 2}
    T_send_eT    :: AbstractArray{Float64, 2}

    eT_MoLap_eT    :: AbstractArray{Float64, 2}

    tool_mtx

    wksp_eT        :: Array
 
    function DiffusionSolver(;
        gi             :: PolelikeCoordinate.GridInfo,
        mask2          :: Union{AbstractArray{Float64, 2}, Nothing} = nothing,
        K              :: Float64,
        Δt             :: Float64,
        M              :: Union{DynamicAdvSpeedUpMatrix, Nothing} = nothing
    )

        α = 1.0 / (K * Δt)

        if M == nothing
            println("Mask not provided. Use mask provided by GridInfo.")
            M = DynamicAdvSpeedUpMatrix(;
                gi = gi,
                Nz = 1,
                mask2 = mask2,
            )
        end

        Nx = gi.Nx
        Ny = gi.Ny
        Nyp1 = Ny+1

        # need a mask excluding bnd points
        # Create coversion matrix and its inverse
        #println(size(M.filter_U))
        #println(Nx, ",", Ny)
        T_num        = reshape(M.filter_T * collect(1:(Nx*Ny)), Nx, Ny)
        active_T_num = T_num[ T_num .!= 0.0 ]
        eT_send_U    = M.op.U_I_U[active_U_num, :]
        U_send_eT    = eT_send_U' |> sparse
        eT_I_eT      = eT_send_U * U_send_eT
        eT_Lap_eT    = eT_send_U * M.U_Lap_U * U_send_eT
        eT_MoLap_eT  = eT_Lap_eT - α * eT_I_eT
        MoLap_eT     = lu(eT_MoLap_eT)

        MoLap_VV     = lu(VV_MoLap_VV)

        tool_mtx = (
            MoLap_VV = MoLap_VV,
            MoLap_eT = MoLap_eT,
        )
        
        eT_length = size(eT_MoLap_eT)[1]
        VV_length = size(VV_MoLap_VV)[1]

        wksp_eT = [ zeros(Float64, eT_length), zeros(Float64, eT_length) ]
        wksp_VV = [ zeros(Float64, VV_length), zeros(Float64, VV_length) ]

        return new(
            eT_length,
            VV_length,
            M,
            K,
            Δt,
            α,
            eT_send_U,
            U_send_eT,
            VV_send_V,
            V_send_VV,
            eT_MoLap_eT,
            VV_MoLap_VV,
            tool_mtx,
            wksp_eT,
            wksp_VV,
        )

    end
end


function solveDiffusion!(
    ds     :: DiffusionSolver,
    grid   :: Symbol,
    input  :: AbstractArray{Float64},
    output :: AbstractArray{Float64},
)

    if grid == :U
        rhs = ds.wksp_eT[1]
        lhs = ds.wksp_eT[2]
        tool = ds.tool_mtx.MoLap_eT
        send    = ds.eT_send_U
        invsend = ds.U_send_eT
    elseif grid == :V
        rhs = ds.wksp_VV[1]
        lhs = ds.wksp_VV[2]
        tool = ds.tool_mtx.MoLap_VV
        send    = ds.VV_send_V
        invsend = ds.V_send_VV
    else
        throw(ErrorException("Unrecognized grid: " * string(grid)))
    end

    mul!(rhs, send, view(input, :))
    rhs .*= -ds.α
    ldiv!(lhs, tool, rhs)
    mul!(view(output, :), invsend, lhs)

end
