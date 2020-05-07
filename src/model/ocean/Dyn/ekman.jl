function assignEkmanFlow!(
    m :: DynModel,
)

    @unpack m

    # Transform input wind stress vector first
    DisplacedPoleCoordinate.project!(
        ev.gi,
        fr.τx_raw,
        fr.τy_raw,
        fr.τx,
        fr.τy,
        direction=:Forward
    )

    H_ek = ev.H[1]
    H_rf = ev.H[2]

    @. co.M̃ = (fr.τx + fr.τy * im) / (ρ_sw * ev.s)
 
    @. co.ṽ_ek =   M̃ / H_ek
    @. co.ṽ_rf = - M̃ / H_rf

    u_ek = getSpace!(co.wksp, :cT)
    v_ek = getSpace!(co.wksp, :cT)
    u_rf = getSpace!(co.wksp, :cT)
    v_rf = getSpace!(co.wksp, :cT)
view(st.u_total, :, :, 1)
    v_ek = view(st.v_total, :, :, 1)
    u_rf = view(st.u_total, :, :, 2)
    v_rf = view(st.v_total, :, :, 2)

    @. u_ek = real(ṽ_ek)
    @. v_ek = imag(ṽ_ek)
 
    @. u_rf = real(ṽ_rf)
    @. v_rf = imag(ṽ_rf)
        
end
