# Load modules used here
using TimerOutputs
using UnicodePlots
using MPI
using CLIMA
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.DGmethods
using CLIMA.DGmethods.NumericalFluxes
using CLIMA.MPIStateArrays
using CLIMA.LowStorageRungeKuttaMethod
using CLIMA.ODESolvers
using CLIMA.GenericCallbacks
using CLIMA.Atmos
using CLIMA.VariableTemplates
using CLIMA.MoistThermodynamics
using CLIMA.AdditiveRungeKuttaMethod
using CLIMA.PlanetParameters
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates
using CLIMA.VTK
using CLIMA.Atmos: vars_state, vars_aux
using DelimitedFiles
using GPUifyLoops
using Random
using CLIMA.IOstrings
using CLIMA.ColumnwiseLUSolver: SingleColumnLU, ManyColumnLU, banded_matrix,
                                banded_matrix_vector_product!
using CLIMA.DGmethods: EveryDirection, HorizontalDirection, VerticalDirection

import ..haspkg 

# Create timeroutput object
to = TimerOutput()

@static if haspkg("CuArrays")
  using CUDAdrv
  using CUDAnative
  using CuArrays
  CuArrays.allowscalar(false)
  const ArrayTypes = (CuArray,)
else
  const ArrayTypes = (Array,)
end

if !@isdefined integration_testing
  const integration_testing =
    parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
end

const seed = MersenneTwister(0)

"""
  Initial condition for Rayleigh Benard problem with Flux boundary conditions
"""
function Initialise_Rising_Bubble!(state::Vars, aux::Vars, (x,y,z), t)
  FT                = eltype(state)
  # can override default gas constants 
  # to moist values later in the driver 
  R_gas::FT         = R_d
  c_p::FT           = cp_d
  c_v::FT           = cv_d
  p0::FT            = MSLP
  gravity::FT       = grav
  # initialise with dry domain 
  q_tot::FT         = 0 
  q_liq::FT         = 0
  q_ice::FT         = 0 
  xc::FT            = 500
  yc::FT            = 500
  zc::FT            = 400
  # perturbation parameters for rising bubble
  r                 = sqrt((x - xc)^2 + (y - yc)^2 + (z-zc)^2)
  rc::FT            = 300
  θ_ref::FT         = 300
  θ_c::FT           = 5.0
  Δθ::FT            = 0.0
  if r <= rc 
    Δθ = θ_c * (1 + cospi(r/rc))/2
  end
  θ                     = θ_ref + Δθ # potential temperature
  π_exner               = 1.0 - gravity / (c_p * θ) * z # exner pressure
  ρ                     = p0 / (R_gas * θ) * (π_exner)^ (c_v / R_gas) # density

  P                     = p0 * (R_gas * (ρ * θ) / p0) ^(c_p/c_v) # pressure (absolute)
  T                     = P / (ρ * R_gas) # temperature
  U, V, W               = 0.0 , 0.0 , 0.0  # momentum components
  # energy definitions
  e_kin                 = (U^2 + V^2 + W^2) / (2*ρ)/ ρ
  e_pot                 = gravity * z
  e_int                 = c_v * (T - T_0) #internal_energy(T, q_tot, q_liq, q_ice)
  E                     = ρ * (e_int + e_kin + e_pot)  #* total_energy(e_kin, e_pot, T, q_tot, q_liq, q_ice)
  state.ρ      = ρ
  state.ρu     = SVector{3,FT}(0,0,0)
  state.ρe     = E
  state.moisture.ρq_tot = FT(0)
end

function run(mpicomm, 
             ArrayType, 
             dim, 
             topl, 
             N, 
             timeend, 
             FT, 
             dt_exp, 
             dt_imex, 
             zmax, 
             zsponge, 
             problem_name, 
             OUTPATH, 
             aspectratio, 
             IsExplicit, 
             LinearModel,
             SolverMethod)

  # Grid setup (topl contains brickrange information)
  grid = DiscontinuousSpectralElementGrid(topl,
                                          FloatType = FT,
                                          DeviceArray = ArrayType,
                                          polynomialorder = N,
                                         )
  # Model definition
  model = AtmosModel(FlatOrientation(),
                     HydrostaticState(IsothermalProfile(FT(300)), FT(0)),
                     SmagorinskyLilly{FT}(0.23),
                     EquilMoist(), 
                     NoRadiation(),
                     Gravity(), 
                     NoFluxBC(),
                     Initialise_Rising_Bubble!)
  # Balancelaw description
  dg = DGModel(model,
               grid,
               Rusanov(),
               CentralNumericalFluxDiffusive(),
               CentralGradPenalty(),
               direction=EveryDirection())

  linmodel = LinearModel(model)
  
  vdg = DGModel(linmodel,
                grid,
                Rusanov(),
                CentralNumericalFluxDiffusive(),
                CentralGradPenalty(),
                auxstate=dg.auxstate,
                direction=VerticalDirection())
  
  Q = init_ode_state(dg, FT(0))
  eng0= norm(Q)

  diagx = [FT(0),]
  diagy = [FT(1),]
  
  # Setup VTK output callbacks
  output_interval = 100
  step = [0]
    cbvtk = GenericCallbacks.EveryXSimulationSteps(output_interval) do (init=false)
    mkpath(OUTPATH)
    outprefix = @sprintf("%s/RB_%dD_mpirank%04d_step%04d", OUTPATH, dim,
                           MPI.Comm_rank(mpicomm), step[1])
    @debug "doing VTK output" outprefix
    writevtk(outprefix, Q, dg, flattenednames(vars_state(model,FT)),
             dg.auxstate, flattenednames(vars_aux(model,FT)))

    step[1] += 1
    nothing
  end

  # ---------------- SOLVER -----------------------#

  if IsExplicit == 1
    numberofsteps = convert(Int64, cld(timeend, dt_exp))
    dt_exp = timeend / numberofsteps
    @info "EXP timestepper" dt_exp numberofsteps dt_exp*numberofsteps timeend
    solver = LSRK54CarpenterKennedy(dg, Q; dt = dt_exp, t0 = 0)
    @timeit to "solve! EX RTB - Explicit LSRK54CarpenterKennedy $aspectratio $dt_exp $timeend" solve!(Q, solver; numberofsteps=numberofsteps, callbacks=(),adjustfinalstep=false)
  else
    numberofsteps = convert(Int64, cld(timeend, dt_imex))
    dt_imex = timeend / numberofsteps
    @info "1DIMEX timestepper" dt_imex numberofsteps dt_imex*numberofsteps timeend
    solver = SolverMethod(dg, vdg, SingleColumnLU(), Q;
                                             dt = dt_imex, t0 = 0,
                                             split_nonlinear_linear=false)
    @timeit to "solve! IMEX RTB- $LinearModel $SolverMethod $aspectratio $dt_imex $timeend" solve!(Q, solver; numberofsteps=numberofsteps, callbacks=(),adjustfinalstep=false)
  end
  
  #Get statistics at the end of the run:
  current_time_str = string(ODESolvers.gettime(solver))
  
  # Print some end of the simulation information
 engf = norm(Q)
  Qe = init_ode_state(dg, FT(timeend))
  engfe = norm(Qe)
  errf = euclidean_distance(Q, Qe)
  @info @sprintf """Finished
  norm(Q) / norm(Q₀)      = %.16e
  -------------------------------
  """ engf/eng0 
  
  return nothing
end

using Test
let
  MPI.Initialized() || MPI.Init()
  mpicomm = MPI.COMM_WORLD
  ll = uppercase(get(ENV, "JULIA_LOG_LEVEL", "INFO"))
  loglevel = ll == "DEBUG" ? Logging.Debug :
    ll == "WARN"  ? Logging.Warn  :
    ll == "ERROR" ? Logging.Error : Logging.Info
  logger_stream = MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull
  global_logger(ConsoleLogger(logger_stream, loglevel))
  @static if haspkg("CUDAnative")
      device!(MPI.Comm_rank(mpicomm) % length(devices()))
  end
 # @testset "$(@__FILE__)" for ArrayType in ArrayTypes
 aspectratios = (1,2,4,8,16,32,50,)
 exp_step = (0,1,)
 LinearModels = (AtmosAcousticGravityLinearModel,)#, AtmosAcousticGravityLinearModel)
 IMEXSolverMethods = (ARK548L2SA2KennedyCarpenter,)
 FloatType = (Float32,)
 for explicit in exp_step
   for ArrayType in ArrayTypes 
     for LinearModel in LinearModels 
       for SolverMethod in IMEXSolverMethods
         for aspectratio in aspectratios
           for FT in FloatType
              # Problem type
              # DG polynomial order
              N = 4
              # User defined domain parameters
              Δh = 35
              Δv = Δh/aspectratio
              xmin, xmax = 0, 1500
              ymin, ymax = 0, 1500
              zmin, zmax = 0, 1500

              grid_resolution = [Δh, Δh, Δv]
              domain_size     = [xmin, xmax, ymin, ymax, zmin, zmax]
              dim = length(grid_resolution)

               brickrange = (grid1d(xmin, xmax, elemsize=FT(grid_resolution[1])*N),
                             grid1d(ymin, ymax, elemsize=FT(grid_resolution[2])*N),
                             grid1d(zmin, zmax, elemsize=FT(grid_resolution[end])*N))
                  
              zsponge = FT(1200.0)

              topl = StackedBrickTopology(mpicomm, brickrange,
                                          periodicity = (true, true, false),
                                          boundary=((0,0),(0,0),(1,2)))

              problem_name = "rising_bubble"
              dt_exp = min(Δv/soundspeed_air(FT(330))/N * FT(0.8), Δh/soundspeed_air(FT(330))/N * FT(0.8))
              dt_imex = Δh/soundspeed_air(FT(330))/N * FT(0.8)
              timeend = FT(1)

              #Create unique output path directory:
              OUTPATH = IOstrings_outpath_name(problem_name, grid_resolution)
              #open diagnostics file and write header:
              mpirank = MPI.Comm_rank(MPI.COMM_WORLD)
              C_smag = FT(0.23)
              if explicit == 1
                @info @sprintf """Starting
                ArrayType                 = %s
                ODE_Solver                = LSRK54CarpenterKennedy
                LinearModel               = None
                dt_exp                    = %.5e
                dt_ratio                  = -N/A-
                Δh/Δv                     = %.5e
                """ ArrayType dt_exp aspectratio
              else
                @info @sprintf """Starting
                ArrayType                 = %s
                ODE_Solver                = %s
                LinearModel               = %s
                dt_exp                    = %.5e
                dt_imp                    = %.5e
                dt_ratio                  = %.3e
                Δh/Δv                     = %.5e
                """ ArrayType SolverMethod LinearModel dt_exp dt_imex dt_imex/dt_exp aspectratio
              end
              result1 = run(mpicomm, 
                            ArrayType, 
                            dim, 
                            topl,
                            N, 
                            timeend, 
                            FT, 
                            dt_exp, 
                            dt_imex, 
                            zmax, 
                            zsponge, 
                            problem_name, 
                            OUTPATH, 
                            aspectratio, 
                            explicit, 
                            LinearModel,
                            SolverMethod)
            end
          end
        end
      end
    end
  end
    show(to, allocations = false, compact = true, sortby=:name)
end

#nothing
