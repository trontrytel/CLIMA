using Test, MPI
import GaussQuadrature
using CLIMA
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Mesh.Geometry
using CLIMA.Mesh.Interpolation
using StaticArrays
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
const seed = MersenneTwister(0)

const ArrayType = CLIMA.array_type()
#------------------------------------------------
if !@isdefined integration_testing
  const integration_testing =
    parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
end
#------------------------------------------------
function run()
  MPI.Initialized() || MPI.Init()

#@testset "LocalGeometry" begin
  FT = Float64
  ArrayType = Array

  xmin, ymin, zmin = 0, 0, 0                   # defining domain extent
  xmax, ymax, zmax = 2000, 400, 2000
  xres, yres, zres = FT(200), FT(200), FT(200) # resolution of interpolation grid

  xgrd = range(xmin, xmax, step=xres) 
  ygrd = range(ymin, ymax, step=yres) 
  zgrd = range(zmin, zmax, step=zres) 

#  Ne        = (20,2,20)
  Ne        = (4,2,4)

  polynomialorder = 4 #8 #4
  #-------------------------
  _x, _y, _z = 12, 13, 14
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
                     Initialise_Test!)

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
  intrp_brck = Interpolation_Brick(grid, xres, yres, zres)
  interpolate_brick!(intrp_brck, Q.data, st_no, polynomialorder)
  #------testing
  Nel = length( grid.topology.realelems )

  error = zeros(FT, Nel) 

  for elno in 1:Nel
    fex = similar(intrp_brck.V[elno])
    fex = fcn( intrp_brck.xg[elno] ./ xmax , intrp_brck.yg[elno] ./ ymax , intrp_brck.zg[elno] ./ zmax )
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
end #function run

#-----taken from Test example
function Initialise_Test!(state::Vars, aux::Vars, (x,y,z), t)
  FT         = eltype(state)
	
	# Dummy variables for initial condition function 
  state.ρ     = FT(0) 
  state.ρu    = SVector{3,FT}(0,0,0)
  state.ρe    = FT(0)
  state.moisture.ρq_tot = FT(0)
end
#------------------------------------------------
run()
