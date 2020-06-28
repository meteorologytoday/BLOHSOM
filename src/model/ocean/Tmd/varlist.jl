function getCoordinateVariable(m::TmdMaster)
        return Dict(
            "f"               => ( m.env.gi.c_f,                                                ("Nx", "Ny",) ),
        )
end

function getCompleteVariableList(m::TmdMaster)

    st = m.state
    
    return Dict(
        "X"               => ( PermutedDimsArray(st.X, (2, 3, 1, 4)),                ("Nx", "Ny", "Nz", "NX") ),
        "X_ML"            => ( st.X_ML,                                              ("Nx", "Ny",       "NX") ),
        "swflx"           => ( st.swflx,                                             ("Nx", "Ny") ),
        "u_U"             => ( PermutedDimsArray(st.u_U, (2, 3, 1)),                 ("Nx", "Ny", "Nz"      ) ),
    )
end

