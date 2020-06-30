
function calCoarseBuoyancyPressure!(
    m :: TmdModel,
)

    bs = m.bridge_state

    if bs == nothing
        return
    end

    calAverage_f2c!(
        bs.va,
        PermutedDimsArray(m.state.B, (2, 3, 1)),
        bs.B_for_dyn,
    )

end

function projVelocity!(
    m :: TmdModel,
)

    bs = m.bridge_state
    projVertical_c2f!(
        bs.va,
        bs.uc_from_dyn,
        PermutedDimsArray(m.state.u_U, (2, 3, 1)),
    )
 
    projVertical_c2f!(
        bs.va,
        bs.vc_from_dyn,
        PermutedDimsArray(m.state.v_V, (2, 3, 1)),
    )
 
end
