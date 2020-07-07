using ArgParse
using JSON
using Formatting
include("BLOHSOM.jl")
include("../../share/constants.jl")

if !isdefined(Main, :REPL)

    function parse_commandline()
        s = ArgParseSettings()
        @add_arg_table s begin

            "--run-days"
                help = "Simulation days."
                arg_type = Int64
                default = 10
     
            "--output-file"
                help = "Output file."
                arg_type = String
                required = true
     
        end

        return parse_args(ARGS, s)
    end



    parsed = parse_commandline()
    print(json(parsed, 4))

else

    run_days = 3
    output_file = "output.nc"

end


run_days = 5
coupling = 12
substeps_dyn = 2
substeps_tmd = 1
rec_frequency = 1


run_days = 2
coupling = 12
substeps_dyn = 2
substeps_tmd = 1
rec_frequency = 12


Δt_day = 86400.0
Δt_coupling = Δt_day / coupling



run_steps = run_days * coupling


#86400.0

Nz_f = 20

hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.ocn.gx3v7.120323.nc"
topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.gx3v7.nc"

#hrgrid_file = "/seley/tienyiah/CESM_domains/test_domains/domain.lnd.fv4x5_gx3v7.091218.nc"
#topo_file = "/seley/tienyiah/CESM_domains/test_domains/topo.fv4x5.nc"

#hrgrid_file = "/seley/tienyiah/CESM_domains/domain.lnd.fv1.9x2.5_gx1v6.090206.nc"

gf = GridFiles.CurvilinearSphericalGridFile(
    hrgrid_file,
    R = Re,
    Ω = Ωe,
)

#=
Ly = 100e3 * 60.0
gf = GridFiles.CylindricalGridFile(;
        R   = Re,
        Ω   = Ωe,
        Nx   = 120,
        Ny   = 60,
        Ly   = Ly,
        lat0 = 0.0 |> deg2rad,
        β    = 2*Ωe / Re,
)
=#
#xcutoff = 1

#ycutoff = 3

#gf.mask                       .= 1
#gf.mask[:, 1:ycutoff]         .= 0 
#gf.mask[:, end-ycutoff+1:end] .= 0 

#gf.mask[1:xcutoff, :]         .= 0 
#gf.mask[end-xcutoff+1:end, :] .= 0 


#gf.mask = 1.0 .- gf.mask
#ycutoff = 2

#gf.mask[:, 1:ycutoff]         .= 0 
#gf.mask[:, end-ycutoff+1:end] .= 0 

gi = PolelikeCoordinate.genGridInfo(gf);


topo = similar(gf.mask)
topo .= -4000
z_bnd_f = collect(Float64, range(0.0, -500.0, length=21))
append!(z_bnd_f, [-575, -800, -1475, -2000, -3000, -4000])
ocn_env = BLOHSOM.OcnEnv(
    hrgrid                = gf,
    topo_file             = topo_file,
    topo_varname          = "topo",
    Δt                    = Δt_coupling,
    substeps_dyn          = substeps_dyn,
    substeps_tmd          = substeps_tmd,
    z_bnd_f               = z_bnd_f,
    #height_level_counts   = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6],
    height_level_counts   = [2, 2, 2, 2, 2, 2, 14],
    NX_passive            = 0,
    deep_threshold        = 4000.0,
    Kh_m_barotropic       = 10000.0,
    Kh_m_baroclinic       = 1000.0,
    Kv_m_baroclinic       = 1e-5,
    τ_barotropic_bottomfric = 120.0 * 86400.0,
    τ_barotropic_coastfric  = 120.0 * 86400.0,
    Kh_X                  = [1e3, 1e-5], 
    Kv_X                  = [1e-5, 1e-5], 
    R                     = 0.48,
    ζ                     = 23.0,
    we_max                = 1e-2,
    MLT_rng               = [10.0, 1000.0],
    X_wr_file             = ["T.nc", "S.nc"],
    X_wr_varname          = ["T", "S"],
    t_X_wr                = [NaN, NaN],
    flow_scheme           = :PROG,
    MLT_scheme            = :prescribe,
    radiation_scheme      = :exponential_decay,
    convective_adjustment = true,
    use_Qflux             = false,
    finding_Qflux         = false;
    mask2                 = gf.mask,
    topo                  = topo,
)

println("##### Initialize model #####")
model = BLOHSOM.init!(ocn_env)

#du[:Φ].data .= 0.01 * exp.(- ( (gi.c_lat * gi.R ).^2 + (gi.R * (gi.c_lon .- π)).^2) / (σ^2.0) / 2)

σz = 150.0
σ = 750e3
σ_warmpool = 1000e3

#du[:X].odata[30:40, 15:25, 1:5, 1] .+= 10.0
#du[:X_ML].odata[30:40, 15:25, 1] .+= 10.0



z_mid = (z_bnd_f[1:end-1] + z_bnd_f[2:end]) / 2
z_mid = repeat(reshape(z_mid, :, 1, 1), outer=(1, gf.Nx, gf.Ny))
mask3 = repeat(reshape(gf.mask, 1, gf.Nx, gf.Ny), outer=(length(z_bnd_f)-1, 1, 1))

basic_T   = z_mid * 0
anomaly_T = z_mid * 0
anomaly_T2 = z_mid * 0

f = (z,)-> (tanh((z-(-200))/50)+1)/2

@. basic_T = 10 + 10 * f(z_mid) #10 + 15 * exp(z_mid / σz)

#warmpool_bnd_z = - 300.0 * exp.(- ( (gi.c_lat * gi.R ).^2 + (gi.R * (gi.c_lon .- π)/2).^2) / (σ_warmpool^2.0) / 2)
#warmpool_bnd_z = - 300.0 * exp.(- ( (gi.c_y .- ( Ly/4.0)).^2 + (gi.R * (gi.c_lon .- π)/2).^2) / (σ_warmpool^2.0) / 2)
#warmpool_bnd_z = - 200.0 * exp.(- ( (gi.c_y *0).^2 + (gi.R * (gi.c_lon .- π)/2).^2) / (σ_warmpool^2.0) / 2)
warmpool_bnd_z = - 200.0 * exp.(- ( (gi.c_lat * gi.R ).^2 + (gi.R * (gi.c_lon .- π)/2).^2) / (σ_warmpool^2.0) / 2)
warmpool_bnd_z = repeat(reshape(warmpool_bnd_z, 1, size(warmpool_bnd_z)...), outer=(size(z_mid)[1], 1, 1))


anomaly_T[z_mid .> warmpool_bnd_z] .= 5.0

#warmpool_bnd_z = - 200.0 * exp.(- ( ((gi.c_lat.- (80* π/180)) * gi.R ).^2 + (gi.R * 0 * (gi.c_lon .- π)/2).^2) / (σ_warmpool^2.0) / 2)
#warmpool_bnd_z = repeat(reshape(warmpool_bnd_z, 1, size(warmpool_bnd_z)...), outer=(size(z_mid)[1], 1, 1))
#anomaly_T2[z_mid .> warmpool_bnd_z] .= 20.0


total_T = basic_T  #+ anomaly_T #+ anomaly_T2
total_T[z_mid .> warmpool_bnd_z] .= 25.0
#total_T[1:4, :, :] .= 20.0

h_ML = gf.mask * 0 .+ 10.0

#total_T .= 10.0
total_T[mask3 .== 0.0] .= 0

tmd_state = model.tmd_engine.state
tmd_state.h_ML .= h_ML
tmd_state.X_ML[:, :, 1] .= total_T[1, :, :]
tmd_state.X[:, :, :, 1] .= total_T

recorder = BLOHSOM.getBasicRecorder(model)
RecordTool.setNewNCFile!(recorder, output_file)
RecordTool.record!(recorder)
RecordTool.avgAndOutput!(recorder)


coupling_cnt = 0
total_cnt = 0
@time for day=1:run_days
    println(format("##### Run Day {:d} #####", day ))

   
#=    if day <= 10 
        du[:SWFLX].data .= - 1000.0 * exp.(- ( (gi.c_lat * gi.R ).^2 + (gi.R * (gi.c_lon .- π)/2).^2) / (σ_warmpool^2.0) / 2)
    else
        du[:SWFLX].data .= 0.0
    end
=#
    @time for c = 1:coupling
        print(format(" - Coupling: {:d}/{:d}\r", c, coupling))
        BLOHSOM.stepModel!(model, false)
        RecordTool.record!(recorder)

        global coupling_cnt += 1
        global total_cnt += 1
        
        if coupling_cnt == rec_frequency
            print("Record!")
            RecordTool.avgAndOutput!(recorder)
            coupling_cnt =  0
        end

    end
        




end
println("Done.")
