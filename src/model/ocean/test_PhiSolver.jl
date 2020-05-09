include("../../share/constants.jl")
include("../../share/PolelikeCoordinate.jl")
include("MatrixOperators.jl")
include("Dyn/AdvectionSpeedUpMatrix_dyn.jl")
include("Dyn/PhiSolver.jl")

H  = 1000.0;
Δt = 86400.0;

n = 6;

println("Create Gridinfo");
Nx = 120
Ny = 60
Ly = 100e3 * 150.0
gf = GridFiles.CylindricalGridFile(;
        R   = Re,
        Ω   = Ωe,
        Nx   = Nx,
        Ny   = Ny,
        Ly   = Ly,
        lat0 = 0.0 |> deg2rad,
        β    = Ωe / Re,
)

gi = PolelikeCoordinate.genGridInfo(gf);


mask2 = ones(Float64, gi.Nx, gi.Ny);

mask2[:, 1] .= 0
mask2[:, end] .= 0

mask2[50:60, 25:35] .= 0

println("Making ΦSolver")
@time cM = PhiSolver(;
    gi    = gi,
    mask2 = mask2,
    α     = 1 / (Δt^2 * H)
);

#F = cM.tool_mtx.MoLap

F = lu(cM.eT_Lap_eT)


f         =   cos.(π/gi.Ly * gi.c_y) .* sin.(2*gi.c_lon);
dfdx_true =   cos.(π/gi.Ly * gi.c_y) .* cos.(2*gi.c_lon) * 2 / gi.R;
dfdy_true = - sin.(π/gi.Ly * gi.c_y) .* sin.(2*gi.c_lon) * π / gi.Ly;
Lapf_true =   - f .* ( (2/ gi.R)^2 + (π/gi.Ly)^2 );

f .*= mask2


f = f[:]

dfdx = cM.M.U_∂x_T * f
dfdy = cM.M.V_∂y_T * f

Lapf = cM.T_Lap_T * f

solve_Lapf = cM.T_send_eT * (F \ (cM.eT_send_T * Lapf))

reshape2 = (m,) -> reshape(m, gi.Nx, :)

f    = reshape2(f)
Lapf = reshape2(Lapf)
solve_Lapf = reshape2(solve_Lapf)
dfdx = reshape2(dfdx)
dfdy = reshape2(dfdy)

using NCDatasets
Dataset("output_phisolver.nc", "c") do ds

    defDim(ds, "Nx", gi.Nx)
    defDim(ds, "Ny", gi.Ny)
    defDim(ds, "Nyp1", gi.Ny+1)

    for (varname, vardata, vardim, attrib) in [
        ("f",  f, ("Nx", "Ny"), Dict()),
        ("dfdx",  dfdx, ("Nx", "Ny"), Dict()),
        ("dfdy",  dfdy, ("Nx", "Nyp1"), Dict()),
        ("Lapf",  Lapf, ("Nx", "Ny"), Dict()),
        ("solve_Lapf",  solve_Lapf, ("Nx", "Ny"), Dict()),
        ("dfdx_true",  dfdx_true, ("Nx", "Ny"), Dict()),
        ("dfdy_true",  dfdy_true, ("Nx", "Ny"), Dict()),
        ("Lapf_true",  Lapf_true, ("Nx", "Ny"), Dict()),
        ("mask2",  mask2, ("Nx", "Ny"), Dict()),

    ]

        if ! haskey(ds, varname)
            var = defVar(ds, varname, Float64, vardim)
            var.attrib["_FillValue"] = 1e20
        end

        var = ds[varname]
        
        for (k, v) in attrib
            var.attrib[k] = v
        end

        rng = []
        for i in 1:length(vardim)-1
            push!(rng, Colon())
        end
        push!(rng, 1:size(vardata)[end])
        var[rng...] = vardata

    end

end

