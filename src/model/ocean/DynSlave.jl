#include("../../share/constants.jl")

mutable struct DynDuke

    model        :: Dyn.DynModel

    ocn_env      :: OcnEnv
    shared_data  :: SharedData

    data_exchanger :: DataExchanger

    function DynDuke(
        ocn_env      :: OcnEnv,
        shared_data  :: SharedData,
    )
       
        gi = PolelikeCoordinate.genGridInfo(
            ocn_env.hrgrid,
        )
       

        z_bnd_c = zeros(Float64, length(ocn_env.height_level_counts) + 1)
        
        z_bnd_c[1] = ocn_env.z_bnd_f[1]
        idx = 1
        for (k, cnt) in enumerate(ocn_env.height_level_counts)
            idx += cnt
            z_bnd_c[k+1] = ocn_env.z_bnd_f[idx]
        end

        if z_bnd_c[end] != ocn_env.z_bnd_f[end]
            throw(ErrorException("Mismatch"))
        end

 
        model = Dyn.DynModel(
            mode    = ocn_env.flow_scheme,
            gi      = gi,
            Δt      = ocn_env.Δt / ocn_env.substeps_dyn,
            Kh_barotropic = ocn_env.Kh_m_barotropic,
            Kh_baroclinic = ocn_env.Kh_m_baroclinic,
            Kv_baroclinic = ocn_env.Kv_m_baroclinic,
            τ_barotropic_bottomfric = ocn_env.τ_barotropic_bottomfric,
            τ_barotropic_coastfric  = ocn_env.τ_barotropic_coastfric,
            z_bnd = z_bnd_c,
            mask    = ocn_env.mask2_deep,
        )         

        data_exchanger = DataExchanger([
            :FR_TMD,
            :TO_TMD,
        ])

        return new(
            model, 
            ocn_env,
            shared_data,
            data_exchanger,
        )

    end

end

