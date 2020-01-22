# Load Packages
using MPI
using CLIMA
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.DGmethods
using CLIMA.DGmethods.NumericalFluxes
using CLIMA.MPIStateArrays
using CLIMA.LowStorageRungeKuttaMethod
using CLIMA.StrongStabilityPreservingRungeKuttaMethod
using CLIMA.MultirateRungeKuttaMethod
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
using CLIMA.DGmethods: EveryDirection, HorizontalDirection, VerticalDirection

const ArrayType = CLIMA.array_type()

if !@isdefined integration_testing
  const integration_testing =
    parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
end

# -------------- Problem constants ------------------- #
const xmin      = 0
const ymin      = 0
const zmin      = 0
const xmax      = 500
const ymax      = 500
const zmax      = 500
const Ne        = (10,10,10)
const polynomialorder = 4
const dim       = 3
const Δx        = (xmax-xmin)/(Ne[1]*polynomialorder+1)
const Δy        = (ymax-ymin)/(Ne[2]*polynomialorder+1)
const Δz        = (zmax-zmin)/(Ne[3]*polynomialorder+1)
const Δ         = cbrt(Δx * Δy * Δz)
const dt        = 0.005
const timeend   = 100
const T_bot     = 299
const T_lapse   = -0.01
const T_top     = T_bot + T_lapse*zmax
const C_smag    = 0.18
# ------------- Initial condition function ----------- #
function initialise_rayleigh_benard!(state::Vars, aux::Vars, (x1,x2,x3), t)
  FT            = eltype(state)
  R_gas::FT     = R_d
  c_p::FT       = cp_d
  c_v::FT       = cv_d
  γ::FT         = c_p / c_v
  p0::FT        = MSLP
  δT            = sinpi(6*x3/(zmax-zmin)) * cospi(6*x3/(zmax-zmin))
  δw            = sinpi(6*x3/(zmax-zmin)) * cospi(6*x3/(zmax-zmin))
  ΔT            = T_lapse * x3 + δT
  T             = T_bot + ΔT
  P             = p0*(T/T_bot)^(grav/R_gas/T_lapse)
  ρ             = P / (R_gas * T)
  ρu, ρv, ρw    = FT(0) , FT(0) , ρ * δw
  E_int         = ρ * c_v * (T-T_0)
  E_pot         = ρ * grav * x3
  E_kin         = ρ * FT(1/2) * δw^2
  ρe_tot        = E_int + E_pot + E_kin
  state.ρ       = ρ
  state.ρu      = SVector(ρu, ρv, ρw)
  state.ρe      = ρe_tot
  state.moisture.ρq_tot = FT(0)
end
# --------------- Driver definition ------------------ #
function run(mpicomm,
             topl, dim, Ne, polynomialorder,
             timeend, FT, dt, model)
  # -------------- Define grid ----------------------------------- #
  grid = DiscontinuousSpectralElementGrid(topl,
                                          FloatType = FT,
                                          DeviceArray = ArrayType,
                                          polynomialorder = polynomialorder
                                           )
  # -------------- Define model ---------------------------------- #
  model = model

  fast_model = AtmosAcousticGravityLinearModel(model)
  
  slow_model = RemainderModel(model,(fast_model,))
  # -------------- Define dgbalancelaw --------------------------- #
  dg = DGModel(model,
               grid,
               Rusanov(),
               CentralNumericalFluxDiffusive(),
               CentralGradPenalty(),
               direction=EveryDirection())

  fast_dg = DGModel(fast_model,
                     grid,
                     Rusanov(),
                     CentralNumericalFluxDiffusive(),
                     CentralGradPenalty();
                     auxstate = dg.auxstate)

  slow_dg = DGModel(slow_model,
                     grid,
                     Rusanov(),
                     CentralNumericalFluxDiffusive(),
                     CentralGradPenalty();
                     auxstate = dg.auxstate)

  Q = init_ode_state(dg, FT(0))

  #solver = LSRK54CarpenterKennedy(dg, Q; dt = dt, t0 = 0)
  dt_fast = FT(0.005)
  dt_slow = FT(0.005)
  #slow_solver = LSRK144NiegemannDiehlBusch(slow_dg, Q; dt = dt_slow, t0 = 0)
  slow_solver = SSPRK33ShuOsher(slow_dg, Q; 
                                dt = dt_slow, t0 = 0)
  fast_solver = SSPRK33ShuOsher(fast_dg, Q; 
                                dt = dt_fast, t0 = 0)
  solver = MultirateRungeKutta((slow_solver,fast_solver))

  eng0 = norm(Q)
  @info @sprintf """Starting
  norm(Q₀) = %.16e
  ArrayType = %s
  FloatType = %s""" eng0 ArrayType FT

  # Set up the information callback (output field dump is via vtk callback: see cbinfo)
  starttime = Ref(now())
  cbinfo = GenericCallbacks.EveryXWallTimeSeconds(60, mpicomm) do (s=false)
    if s
      starttime[] = now()
    else
      energy = norm(Q)
      @info @sprintf("""Update
                     simtime = %.16e
                     runtime = %s
                     norm(Q) = %.16e""", ODESolvers.gettime(solver),
                     Dates.format(convert(Dates.DateTime,
                                          Dates.now()-starttime[]),
                                  Dates.dateformat"HH:MM:SS"),
                     energy)
    end
  end

  step = [0]
  cbvtk = GenericCallbacks.EveryXSimulationSteps(3000)  do (init=false)
    mkpath("./vtk-rb/")
      outprefix = @sprintf("./vtk-rb/RB_%dD_mpirank%04d_step%04d", dim,
                           MPI.Comm_rank(mpicomm), step[1])
      @debug "doing VTK output" outprefix
      writevtk(outprefix, Q, dg, flattenednames(vars_state(model,FT)), dg.auxstate, flattenednames(vars_aux(model,FT)))
      step[1] += 1
      nothing
  end

  solve!(Q, solver; timeend=timeend, callbacks=(cbinfo,cbvtk))
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

  for FT in (Float32,)
    SGSmodels = (AnisoMinDiss{FT}(1), Vreman{FT}(C_smag), SmagorinskyLilly{FT}(C_smag))
    Expected = (FT(9.9859344959259033e-01),FT(1.0038942098617554e+00),FT(1.0027571916580200e+00))
    T_min = FT(289)
    T_s = FT(290.4)
    Γ_lapse = FT(grav/cp_d)
    Temp = LinearTemperatureProfile(T_min,T_s,Γ_lapse)
    RelHum = FT(0)
    for ii=1:length(SGSmodels)
      model = AtmosModel(FlatOrientation(),
                         HydrostaticState(Temp,RelHum),
                         SGSmodels[ii],
                         EquilMoist(),
                         NoRadiation(),
                         NoSubsidence{FT}(),
                         Gravity(),
                         RayleighBenardBC{FT}(T_bot,T_top),
                         initialise_rayleigh_benard!)
      brickrange = (range(FT(xmin); length=Ne[1]+1, stop=xmax),
                    range(FT(ymin); length=Ne[2]+1, stop=ymax),
                    range(FT(zmin); length=Ne[3]+1, stop=zmax))
      topl = StackedBrickTopology(mpicomm, brickrange, periodicity = (true, true, false), boundary=((0,0),(0,0),(1,2)))
      engf_eng0 = run(mpicomm,
                      topl, dim, Ne, polynomialorder,
                      timeend, FT, dt, model)
      @test engf_eng0 ≈ Expected[ii]
    end
  end
end

#nothing
