mutable struct BridgeState

    va          :: VerticalAverager
    B_to_dyn    :: Union{AbstractArray{Float64, 3}, Nothing}
    uc_from_dyn :: Union{AbstractArray{Float64, 3}, Nothing}
    vc_from_dyn :: Union{AbstractArray{Float64, 3}, Nothing}
    
    function BridgeState(
        z_bnd_f :: AbstractArray{Float64, 1},           # z boundaries of fine grid
        height_level_counts   :: AbstractArray{Int64, 1},
    )
        return new(
            VerticalAverager(z_bnd_f=z_bnd_f, height_level_counts=height_level_counts), 
            nothing,
            nothing,
            nothing,
        )
    end

end



