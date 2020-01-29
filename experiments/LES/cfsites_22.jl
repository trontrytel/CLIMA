using Distributions
using Random
using StaticArrays
using Test
using Printf
using NCDatasets
using Dierckx
# ------------------------------------
using CLIMA
using CLIMA.Atmos
using CLIMA.GenericCallbacks
using CLIMA.LowStorageRungeKuttaMethod
using CLIMA.Mesh.Filters
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters
using CLIMA.VariableTemplates
# ------------------------------------
const ArrayType = CLIMA.array_type()
# ------------------------------------
if !@isdefined integration_testing
  const integration_testing =
    parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
end


# ------------------------------------
# ---------------- Get initial condition from NCData ------------------------------ # 
function get_ncdata()
  data = Dataset("/home/asridhar/CLIMA/datasets/cfsites_forcing.2010071518.nc","r");
  # Load specific site group via numeric ID in NetCDF file (requires generalisation)
  siteid = data.group["site22"];
  # Allow strings to be read as varnames
  function str2var(str::String, var::Any)
    str = Symbol(str);
    @eval(($str)=($var));
  end
  # Load all variables
  for (varname,var) in siteid
    str2var(varname,var[:,1]);
  end
  initdata = [height pfull temp ucomp vcomp sphum]
  return initdata
end

# ------------- Initial condition function ----------- #
"""
CMIP6 Test Dataset - cfsites
@Article{gmd-10-359-2017,
AUTHOR = {Webb, M. J. and Andrews, T. and Bodas-Salcedo, A. and Bony, S. and Bretherton, C. S. and Chadwick, R. and Chepfer, H. and Douville, H. and Good, P. and Kay, J. E. and Klein, S. A. and Marchand, R. and Medeiros, B. and Siebesma, A. P. and Skinner, C. B. and Stevens, B. and Tselioudis, G. and Tsushima, Y. and Watanabe, M.},
TITLE = {The Cloud Feedback Model Intercomparison Project (CFMIP) contribution to CMIP6},
JOURNAL = {Geoscientific Model Development},
VOLUME = {10},
YEAR = {2017},
NUMBER = {1},
PAGES = {359--384},
URL = {https://www.geosci-model-dev.net/10/359/2017/},
DOI = {10.5194/gmd-10-359-2017}
}
"""
function init_cfsites!(state::Vars, 
                            aux::Vars, 
                            (x1,x2,x3), 
                            t,
                            extra_args)

  # Interpolate from CFSite Data
  FT = eltype(state)

  (spl_temp, spl_pfull, spl_ucomp, spl_vcomp, spl_sphum) = extra_args

  T     = FT(spl_temp(x3))
  q_tot = FT(spl_sphum(x3))
  u     = FT(spl_ucomp(x3))
  v     = FT(spl_vcomp(x3))
  P     = FT(spl_pfull(x3))

  ρ     = air_density(T,P,PhasePartition(q_tot))
  e_int = internal_energy(T,PhasePartition(q_tot))
  e_kin = (u^2 + v^2)/2  
  e_pot = grav * x3 
  # Assignment of state variables
  state.ρ = ρ
  state.ρu = ρ * SVector(u,v,0)
  state.ρe = ρ * (e_kin + e_pot + e_int)
  state.moisture.ρq_tot = ρ * q_tot
end

function config_cfsites(FT, N, resolution, xmax, ymax, zmax)
  # Moisture
  moisture = EquilMoist()
  # Radiation model
  radiation = NoRadiation()
  # Subsidence model
  subsidence = NoSubsidence{FT}()
  # Source
  source = Gravity()
  # Sponge
  c_sponge = 1
  zsponge = FT(zmax * 0.80)
  u_relaxation = SVector(0,0,0)
  rayleigh_sponge = RayleighSponge{FT}(zmax, zsponge, c_sponge, u_relaxation, 4)
  # SGS Filter Constant
  C_smag = FT(0.21)
  # Boundary Conditions
  bc = NoFluxBC()

  config = CLIMA.LES_Configuration("CFSite_Demo", N, resolution, xmax, ymax, zmax,
                                   init_cfsites!,
                                   solver_type=CLIMA.ExplicitSolverType(LSRK144NiegemannDiehlBusch),
                                   C_smag=C_smag,
                                   moisture=moisture,
                                   radiation=radiation,
                                   subsidence=subsidence,
                                   sources=source,
                                   bc=bc)
  return config
end

function main()
    CLIMA.init()
    FT = Float64
    # DG polynomial order
    N = 4
    # Domain resolution and size
    Δh = FT(200)
    Δv = FT(100)
    resolution = (Δh, Δh, Δv)
    xmax = 9600
    ymax = 9600
    zmax = 6000
    t0 = FT(0)
    timeend = FT(0.2)

    # CFSite Specific Information
    initdata = get_ncdata()
    
    z = initdata[:,1];
    pfull = initdata[:,2];
    temp  = initdata[:,3];
    ucomp = initdata[:,4];
    vcomp = initdata[:,5];
    sphum = initdata[:,6];

    splines = (spl_temp = Spline1D(z,temp), 
               spl_pfull=Spline1D(z,pfull), 
               spl_ucomp=Spline1D(z,ucomp),
               spl_vcomp=Spline1D(z,vcomp),
               spl_sphum=Spline1D(z,sphum))

    driver_config = config_cfsites(FT, 
                                   N, resolution, 
                                   xmax, ymax, zmax)

    solver_config = CLIMA.setup_solver(t0, timeend, 
                                       driver_config, 
                                       extra_args=splines; 
                                       forcecpu=true)
    
    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(2) do (init=false)
        Filters.apply!(solver_config.Q, 6, solver_config.dg.grid, TMARFilter())
        nothing
    end
    result = CLIMA.invoke!(solver_config;
                          user_callbacks=(cbtmarfilter,cbinformation),
                          check_euclidean_distance=true)
end
main()
