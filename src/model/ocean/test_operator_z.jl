include("../../share/constants.jl")
#include("../../share/GridFiles.jl")
include("../../share/PolelikeCoordinate.jl")

include("MatrixOperators.jl")
include("Tmd/AdvectionSpeedUpMatrix.jl")



gf = GridFiles.CylindricalGridFile(;
        R   = Re,
        Ω   = Ωe,
        Nx   = 5,
        Ny   = 3,
        Ly   = 100e3 * 60,
        lat0 = 0.0 |> deg2rad,
        β    = Ωe / Re,
)

gi = PolelikeCoordinate.genGridInfo(gf);

gf.mask .= 1
Nz = 10


Nz_av = zeros(Int64, gi.Nx, gi.Ny)
Nz_av .= Nz

mask3 = ones(Float64, Nz, gi.Nx, gi.Ny)
noflux_x_mask3 = ones(Float64, Nz, gi.Nx, gi.Ny)
noflux_y_mask3 = ones(Float64, Nz, gi.Nx, gi.Ny+1)

Δz_T = zeros(Float64, Nz, gi.Nx, gi.Ny)
Δz_T .= 10.0


Δz_W = zeros(Float64, Nz+1, gi.Nx, gi.Ny)
Δz_W .= 10


println("Making DynamicAdvSpeedUpMatrix")
M = AdvectionSpeedUpMatrix(;
    gi       = gi,
    Nz       = Nz,
    Nz_av    = Nz_av,
    mask3    = mask3,
    noflux_x_mask3 = noflux_x_mask3,
    noflux_y_mask3 = noflux_y_mask3,
    Δz_T = Δz_T,
    Δz_W = Δz_W,
);

f_W = zeros(Nz+1, gi.Nx, gi.Ny)
f_W[:] = collect(1:length(f_W))


f_W_filtered = reshape( M.filter_W * f_W[:], size(f_W)... )





