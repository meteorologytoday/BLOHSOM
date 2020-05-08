function getCoordinateVariable(m::DynModel)
        return Dict(
            "f"               => ( m.env.gi.c_f,                         ("Nx", "Ny",) ),
        )
end

function getCompleteVariableList(m::DynModel)
        s = m.state
        c = m.core
        f = m.forcing
        return Dict(
            "Phi"               => ( s.Φ,                                ("Nx", "Ny",) ),
            "U"                 => ( s.U,                                ("Nx", "Ny",) ),
            "V"                 => ( s.V,                                ("Nx", "Ny",) ),
            "u_total"           => ( s.u_total,                          ("Nx", "Ny", "Nz") ),
            "v_total"           => ( s.v_total,                          ("Nx", "Ny", "Nz") ),
            "u"                 => ( s.u,                                ("Nx", "Ny", "Nz") ),
            "v"                 => ( s.v,                                ("Nx", "Ny", "Nz") ),
            "B"                 => ( f.B,                                ("Nx", "Ny", "Nz") ),
        )
end

