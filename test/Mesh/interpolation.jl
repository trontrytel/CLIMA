using Test, MPI
import GaussQuadrature
using CLIMA
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Mesh.Geometry
using CLIMA.Mesh.Interpolation
using StaticArrays

#using Plots
#------------------------------------------------
using CLIMA.DGmethods
using CLIMA.DGmethods.NumericalFluxes
using CLIMA.MPIStateArrays
using CLIMA.LowStorageRungeKuttaMethod
using CLIMA.ODESolvers
using CLIMA.GenericCallbacks
using CLIMA.Atmos
using CLIMA.VariableTemplates
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters
using CLIMA.TicToc
using LinearAlgebra
using StaticArrays
using Logging, Printf, Dates
using CLIMA.VTK

using CLIMA.Atmos: vars_state, vars_aux

using Random
using Statistics
const seed = MersenneTwister(0)

const ArrayType = CLIMA.array_type()


#------------------------------------------------
if !@isdefined integration_testing
  const integration_testing =
    parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
end
#------------------------------------------------
function run_brick_interpolation_test()
  MPI.Initialized() || MPI.Init()

#@testset "LocalGeometry" begin
  FT = Float64
  ArrayType = Array

  xmin, ymin, zmin = 0, 0, 0                   # defining domain extent
  xmax, ymax, zmax = 2000, 400, 2000
  xres = [FT(200), FT(200), FT(200)] # resolution of interpolation grid

  xgrd = range(xmin, xmax, step=xres[1]) 
  ygrd = range(ymin, ymax, step=xres[2]) 
  zgrd = range(zmin, zmax, step=xres[3]) 

#  Ne        = (20,2,20)
  Ne        = (4,2,4)

  polynomialorder = 8 #8 #4
  #-------------------------
  _x, _y, _z = CLIMA.Mesh.Grids.vgeoid.x1id, CLIMA.Mesh.Grids.vgeoid.x2id, CLIMA.Mesh.Grids.vgeoid.x3id
  _ρ, _ρu, _ρv, _ρw = 1, 2, 3, 4
  #-------------------------

  brickrange = (range(FT(xmin); length=Ne[1]+1, stop=xmax),
                range(FT(ymin); length=Ne[2]+1, stop=ymax),
                range(FT(zmin); length=Ne[3]+1, stop=zmax))
  topl = StackedBrickTopology(MPI.COMM_SELF, brickrange, periodicity = (true, true, false))

  grid = DiscontinuousSpectralElementGrid(topl,
                                          FloatType = FT,
                                          DeviceArray = ArrayType,
                                          polynomialorder = polynomialorder
                                          )

  model = AtmosModel(FlatOrientation(),
                     NoReferenceState(),
					 ConstantViscosityWithDivergence(FT(0)),
                     EquilMoist(),
                     NoRadiation(),
                     (Gravity()),
					 NoFluxBC(),
                     Initialize_Brick_Interpolation_Test!)

  dg = DGModel(model,
               grid,
               Rusanov(),
               CentralNumericalFluxDiffusive(),
               CentralGradPenalty())

  Q = init_ode_state(dg, FT(0))
  #------------------------------
  x1 = @view grid.vgeo[:,_x,:]
  x2 = @view grid.vgeo[:,_y,:]
  x3 = @view grid.vgeo[:,_z,:]
 
  st_no = _ρ # state vector
  elno = 10

  var = @view Q.data[:,st_no,:]

  #fcn(x,y,z) = x .* y .* z # sample function
  fcn(x,y,z) = sin.(x) .* cos.(y) .* cos.(z) # sample function

  var .= fcn( x1 ./ xmax, x2 ./ ymax, x3 ./ zmax )
  #----calling interpolation function on state variable # st_no--------------------------
  intrp_brck = Interpolation_Brick(grid, xres)
  interpolate_brick_mf!(intrp_brck, Q.data, st_no, polynomialorder)
  #------testing
  Nel = length( grid.topology.realelems )

  error = zeros(FT, Nel) 

  for elno in 1:Nel
    fex = similar(intrp_brck.V[elno])
    fex = fcn( intrp_brck.x[elno][1,:] ./ xmax , intrp_brck.x[elno][2,:] ./ ymax , intrp_brck.x[elno][3,:] ./ zmax )
    error[elno] = maximum(abs.(intrp_brck.V[elno][:]-fex[:]))
  end

  println("==============================================")
  println("l_infinity interpolation error in each element")
  display(error)
  l_infinity_domain = maximum(error)
  println("l_infinity interpolation error in domain")
  display(l_infinity_domain)
  println("==============================================")
  #----------------
end #function run_brick_interpolation_test

#-----taken from Test example
function Initialize_Brick_Interpolation_Test!(state::Vars, aux::Vars, (x,y,z), t)
  FT         = eltype(state)
	
  # Dummy variables for initial condition function 
  state.ρ     = FT(0) 
  state.ρu    = SVector{3,FT}(0,0,0)
  state.ρe    = FT(0)
  state.moisture.ρq_tot = FT(0)
end
#------------------------------------------------
# Cubed sphere, lat/long interpolation test
#----------------------------------------------------------------------------
function run_cubed_sphere()
  CLIMA.init()

  FT = Float64
  mpicomm = MPI.COMM_WORLD

  ll = uppercase(get(ENV, "JULIA_LOG_LEVEL", "INFO"))
  loglevel = Dict("DEBUG" => Logging.Debug,
                  "WARN"  => Logging.Warn,
                  "ERROR" => Logging.Error,
                  "INFO"  => Logging.Info)[ll]

  logger_stream = MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull
  global_logger(ConsoleLogger(logger_stream, loglevel))

  polynomialorder = 4#1#4 #5
  numelem_horz = 3#4 #6
  numelem_vert = 4#1 #1 #1#6 #8

  #-------------------------
  _x, _y, _z = CLIMA.Mesh.Grids.vgeoid.x1id, CLIMA.Mesh.Grids.vgeoid.x2id, CLIMA.Mesh.Grids.vgeoid.x3id
  _ρ, _ρu, _ρv, _ρw = 1, 2, 3, 4
  #-------------------------
#  vert_range = grid1d(FT(planet_radius), FT(planet_radius + setup.domain_height), nelem = numelem_vert)
  vert_range = grid1d(FT(1.0), FT(2.0), nelem = numelem_vert)
  nvert = 4
  
  lat_res  = 5 * π / 180.0 # 5 degree resolution
  long_res = 5 * π / 180.0 # 5 degree resolution
  r_res    = (vert_range[end] - vert_range[1])/FT(nvert) #1000.00    # 1000 m vertical resolution

  #----------------------------------------------------------
  setup = TestSphereSetup{FT}()

  topology = StackedCubedSphereTopology(mpicomm, numelem_horz, vert_range)

  grid = DiscontinuousSpectralElementGrid(topology,
                                          FloatType = FT,
                                          DeviceArray = ArrayType,
                                          polynomialorder = polynomialorder,
                                          meshwarp = CLIMA.Mesh.Topologies.cubedshellwarp)

  model = AtmosModel(SphericalOrientation(),
                     NoReferenceState(),
                     ConstantViscosityWithDivergence(FT(0)),
                     DryModel(),
                     NoRadiation(),
                     nothing, 
                     NoFluxBC(),
                     setup)

  dg = DGModel(model, grid, Rusanov(),
               CentralNumericalFluxDiffusive(), CentralGradPenalty())

  Q = init_ode_state(dg, FT(0))
  #------------------------------
  x1 = @view grid.vgeo[:,_x,:]
  x2 = @view grid.vgeo[:,_y,:]
  x3 = @view grid.vgeo[:,_z,:]
  #------------------------------
  x1_un = @view grid.topology.elemtocoord[1,:,:] # elemtocoord[x/y/z,,vert#, elem#] 
  x2_un = @view grid.topology.elemtocoord[2,:,:] 
  x3_un = @view grid.topology.elemtocoord[3,:,:]


  intrp_cs = Interpolation_Cubed_Sphere(grid, collect(vert_range), lat_res, long_res, r_res)

  elno = 1

  X1 = x1_un[:,elno];  X2 = x2_un[:,elno];  X3 = x3_un[:,elno]; # vertices of el # elno

  xp = [-1.120, 0.7, -0.6]

  xx1 = vcat(X1,xp[1]); xx2 = vcat(X2,xp[2]); xx3 = vcat(X3,xp[3]) 

  println("xp1 = ", xp[1], "; xp2 = ", xp[2], "; xp3 = ", xp[3])

  ξ = invert_trilear_mapping_hex(X1, X2, X3, xp)

  println("ξ = "); display(ξ)
  #----------------------------------------------------------
#  scatter( x1_un[:], x2_un[:], x3_un[:], legend = false)

#  scatter( intrp_cs.x1_uw_grd[:], intrp_cs.x2_uw_grd[:], intrp_cs.x3_uw_grd[:], xlabel = "X", ylabel = "Y", zlabel = "Z", legend = false)

#  xx = intrp_cs.x1_uw_grd[:]
#  yy = intrp_cs.x2_uw_grd[:]
#  zz = intrp_cs.x3_uw_grd[:]

#toler = 1e-12
#loc = (abs.(xx .- 1) .< toler) .* (abs.(xx .+ 1) .> toler) .* 
#      (abs.(yy .- 1) .< toler) .* (abs.(yy .+ 1) .> toler) .* 
#      (abs.(zz .- 1) .< toler) .* (abs.(zz .+ 1) .> toler);

#=  for k in 1:numelem_horz
    println("=======k = ",k,"====================================")
    for j in 1:numelem_vert
    i = j + (k-1)*numelem_vert
    println("xm = ", sum(x1_un[:,i])/length(x1_un[:,i]), "; ym = ", sum(x2_un[:,i])/length(x2_un[:,i]),  "; zm = ", sum(x3_un[:,i])/length(x3_un[:,i] ) )
    end
    println("===========================================")
  end
=#
#  for i in 1:length(grid.topology.realelems)
#    println(i,"). xm = ",mean(x1[:,i]), "; ym = ", mean(x2[:,i]), "; zm = ",  mean(x3[:,i]))
#    if mod(i, length(grid.topology.realelems)/6  ) == 0
#        println("----------------------------------------------------")
#    end
    
#  end
#----------------------------------------------------------------------------

end 
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
Base.@kwdef struct TestSphereSetup{FT}
  p_ground::FT = MSLP
  T_initial::FT = 255
  domain_height::FT = 30e3
end

#----------------------------------------------------------------------------
function (setup::TestSphereSetup)(state, aux, coords, t) 
  # callable to set initial conditions
  FT = eltype(state)

  r = norm(coords, 2)
  h = r - FT(planet_radius)

  scale_height = R_d * setup.T_initial / grav
  p = setup.p_ground * exp(-h / scale_height)

  state.ρ = air_density(setup.T_initial, p)
  state.ρu = SVector{3, FT}(0, 0, 0)
  state.ρe = state.ρ * (internal_energy(setup.T_initial) + aux.orientation.Φ)
  nothing
end
#----------------------------------------------------------------------------


run_brick_interpolation_test()
#run_cubed_sphere()
#------------------------------------------------

