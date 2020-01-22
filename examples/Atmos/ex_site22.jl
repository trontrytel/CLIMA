# Load Packages
using MPI
using CLIMA
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Mesh.Geometry
using CLIMA.DGmethods
using CLIMA.DGmethods.NumericalFluxes
using CLIMA.MPIStateArrays
using CLIMA.LowStorageRungeKuttaMethod
using CLIMA.SubgridScaleParameters
using CLIMA.ODESolvers
using CLIMA.GenericCallbacks
using CLIMA.Atmos
using CLIMA.VariableTemplates
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates
using CLIMA.VTK
using Random
using CLIMA.Atmos: vars_state, vars_aux

using NCDatasets
using Dierckx

const ArrayType = CLIMA.array_type()

if !@isdefined integration_testing
  const integration_testing =
    parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
end

# ------------ SITE 22 Domain Description ----------- # 
const (xmin,xmax)      = (0,9600)
const (ymin,ymax)      = (0,9600)
const (zmin,zmax)      = (0,6000)

const Ne        = (10,10,25)
const polynomialorder = 4

const dim       = 3

const dt        = 0.04
const timeend   = dt

# ---------------- Get initial condition from NCData ------------------------------ # 
function get_ncdata()
  data = Dataset("./cfsites_forcing.2010071518.nc","r");

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

# -------------------------------------------------------------------------------- #

# ------------- Initial condition function ----------- #
"""
CMIP6 Test Dataset - Site22
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

initdata = get_ncdata()

# --------------- Driver definition ------------------ #
function run(mpicomm,
             topl, dim, Ne, polynomialorder,
             timeend, FT, dt)
  # -------------- Define grid ----------------------------------- #
  grid = DiscontinuousSpectralElementGrid(topl,
                                          FloatType = FT,
                                          DeviceArray = ArrayType,
                                          polynomialorder = polynomialorder
                                           )
  # -------------- Define model ---------------------------------- #
  function Initialise_Site22!(state::Vars, aux::Vars, (x1,x2,x3), t,
                              spl_pfull,
                              spl_temp,
                              spl_ucomp, 
                              spl_vcomp,
                              spl_sphum)
    
    FT = eltype(state)
    
    T     = FT(spl_temp(x3))
    q_tot = FT(spl_sphum(x3))
    u     = FT(spl_ucomp(x3))
    v     = FT(spl_vcomp(x3))
    P     = FT(spl_pfull(x3))
    ρ     = air_density(T,P,PhasePartition(q_tot))
    e_int = internal_energy(T,PhasePartition(q_tot))
    e_kin = (u^2 + v^2)/2  
    e_pot = grav * x3 
    
    state.ρ = ρ
    state.ρu = ρ * SVector(u,v,0)
    state.ρe = ρ * (e_kin + e_pot + e_int)
    state.moisture.ρq_tot = ρ * q_tot

  end

  # Get NetCDF data for site 22 and load into initial condition prescription
  initdata = get_ncdata()

  z         = initdata[:, 1];
  pfull     = initdata[:, 2];
  temp      = initdata[:, 3];
  ucomp     = initdata[:, 4];
  vcomp     = initdata[:, 5];
  sphum     = initdata[:, 6];
  
  spl_pfull = Spline1D(z, pfull;k=1)
  spl_temp  = Spline1D(z, temp;k=1)
  spl_ucomp = Spline1D(z, ucomp;k=1)
  spl_vcomp = Spline1D(z, vcomp;k=1)
  spl_sphum = Spline1D(z, sphum;k=1)
  
  Initialise_Site22!(state, aux,(x1,x2,x3),t...) = Initialise_Site22!(state::Vars, aux::Vars, (x1,x2,x3), t,
                                                                                 spl_pfull,
                                                                                 spl_temp,
                                                                                 spl_ucomp, 
                                                                                 spl_vcomp,
                                                                                 spl_sphum)

  model = AtmosModel(FlatOrientation(),
                     NoReferenceState(),
                     SmagorinskyLilly{FT}(0.23),
                     EquilMoist(),
                     NoRadiation(),
                     NoSubsidence{FT}(),
                     Gravity(),
                     NoFluxBC(),
                     Initialise_Site22!)
  
  # -------------- Define dgbalancelaw --------------------------- #
  dg = DGModel(model,
               grid,
               Rusanov(),
               CentralNumericalFluxDiffusive(),
               CentralGradPenalty())


  Q = init_ode_state(dg, FT(0); forcecpu=true)

  lsrk = LSRK54CarpenterKennedy(dg, Q; dt = dt, t0 = 0)

  eng0 = norm(Q)
  @info @sprintf """Starting
  norm(Q₀) = %.16e
  ArrayType = %s
  FloatType = %s""" eng0 ArrayType FT

  # Set up the information callback (output field dump is via vtk callback: see cbinfo)
  starttime = Ref(now())
  cbinfo = GenericCallbacks.EveryXWallTimeSeconds(10, mpicomm) do (s=false)
    if s
      starttime[] = now()
    else
      energy = norm(Q)
      @info @sprintf("""Update
                     simtime = %.16e
                     runtime = %s
                     norm(Q) = %.16e""", ODESolvers.gettime(lsrk),
                     Dates.format(convert(Dates.DateTime,
                                          Dates.now()-starttime[]),
                                  Dates.dateformat"HH:MM:SS"),
                     energy)
    end
  end

  step = [0]
  cbvtk = GenericCallbacks.EveryXSimulationSteps(3000)  do (init=false)
    mkpath("./vtk-cmip6_site22/")
      outprefix = @sprintf("./vtk-cmip6_site22/site22_%dD_mpirank%04d_step%04d", dim,
                           MPI.Comm_rank(mpicomm), step[1])
      @debug "doing VTK output" outprefix
      writevtk(outprefix, Q, dg, flattenednames(vars_state(model,FT)), dg.auxstate, flattenednames(vars_aux(model,FT)))
      step[1] += 1
      nothing
  end

  solve!(Q, lsrk; timeend=timeend, callbacks=(cbinfo,cbvtk))
  # End of the simulation information
  engf = norm(Q)
  Qe = init_ode_state(dg, FT(timeend))
  engfe = norm(Qe)
  errf = euclidean_distance(Q, Qe)
  @info @sprintf """Finished
  norm(Q)                 = %.16e
  norm(Q) / norm(Q₀)      = %.16e
  norm(Q) - norm(Q₀)      = %.16e
  norm(Q - Qe)            = %.16e
  norm(Q - Qe) / norm(Qe) = %.16e
  """ engf engf/eng0 engf-eng0 errf errf / engfe
engf/eng0
end
# --------------- Test block / Loggers ------------------ #
using Test
let
  CLIMA.init()
  mpicomm = MPI.COMM_WORLD
  ll = uppercase(get(ENV, "JULIA_LOG_LEVEL", "INFO"))
  loglevel = ll == "DEBUG" ? Logging.Debug :
    ll == "WARN"  ? Logging.Warn  :
    ll == "ERROR" ? Logging.Error : Logging.Info
  logger_stream = MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull
  global_logger(ConsoleLogger(logger_stream, loglevel))
  FT = Float32
  brickrange = (range(FT(xmin); length=Ne[1]+1, stop=xmax),
                range(FT(ymin); length=Ne[2]+1, stop=ymax),
                range(FT(zmin); length=Ne[3]+1, stop=zmax))
  topl = StackedBrickTopology(mpicomm, brickrange, periodicity = (false, true, false))
  engf_eng0 = run(mpicomm,
                  topl, dim, Ne, polynomialorder,
                  timeend, FT, dt)
end

#nothing
