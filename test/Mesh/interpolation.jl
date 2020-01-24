using Test, MPI
import GaussQuadrature
using CLIMA
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Mesh.Geometry
using CLIMA.Mesh.Interpolation
using StaticArrays
using GPUifyLoops

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
#if !@isdefined integration_testing
#    const integration_testing =
#    parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
#end
#------------------------------------------------
function run_brick_interpolation_test()
  MPI.Initialized() || MPI.Init()

#@testset "LocalGeometry" begin
    FT = Float64
#    ArrayType = ArrayType #Array
    mpicomm = MPI.COMM_WORLD

    xmin, ymin, zmin = 0, 0, 0                   # defining domain extent
    xmax, ymax, zmax = 2000, 400, 2000
    xres = [FT(200), FT(200), FT(200)] # resolution of interpolation grid

#    Ne        = (4,4,4)
#    Ne        = (20,2,20)
#    Ne        = (60,8,60)
    Ne        = (4,2,4)

    polynomialorder = 5 #8# 5 #8 #4
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
                                            polynomialorder = polynomialorder)

    model = AtmosModel(FlatOrientation(),
                     NoReferenceState(),
					 ConstantViscosityWithDivergence(FT(0)),
                     EquilMoist(),
                     NoRadiation(),
                     NoSubsidence{FT}(),
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
 
    st_idx = _ρ # state vector

    var = @view Q.data[:,st_idx,:]

    #fcn(x,y,z) = x .* y .* z # sample function
    fcn(x,y,z) = sin.(x) .* cos.(y) .* cos.(z) # sample function

    var .= fcn( x1 ./ xmax, x2 ./ ymax, x3 ./ zmax )
    #----calling interpolation function on state variable # st_idx--------------------------
    #intrp_brck = InterpolationBrick(grid, xres, FT)
    intrp_brck = InterpolationBrick(grid, xres)


    interpolate_brick!(intrp_brck.offset, intrp_brck.m1_r, intrp_brck.wb, 
                       intrp_brck.ξ1, intrp_brck.ξ2, intrp_brck.ξ3, 
                       intrp_brck.o1, intrp_brck.o2, intrp_brck.o3, 
                       intrp_brck.flg1, intrp_brck.flg2, intrp_brck.flg3, 
                       intrp_brck.fac1, intrp_brck.fac2, intrp_brck.fac3, 
                       intrp_brck.v, Q.data, st_idx)

    for i in 1:20
        @time interpolate_brick!(intrp_brck.offset, intrp_brck.m1_r, intrp_brck.wb, 
                               intrp_brck.ξ1, intrp_brck.ξ2, intrp_brck.ξ3, 
                               intrp_brck.o1, intrp_brck.o2, intrp_brck.o3, 
                               intrp_brck.flg1, intrp_brck.flg2, intrp_brck.flg3, 
                               intrp_brck.fac1, intrp_brck.fac2, intrp_brck.fac3, 
                               intrp_brck.v, Q.data, st_idx)
    end
    #------testing
    Nel = length( grid.topology.realelems )

    error = zeros(FT, Nel) 
    for elno in 1:Nel
        st  = intrp_brck.offset[elno] + 1
        en  = intrp_brck.offset[elno+1]
        if en ≥ st
            fex = fcn( intrp_brck.x[elno][1,:] ./ xmax , intrp_brck.x[elno][2,:] ./ ymax , intrp_brck.x[elno][3,:] ./ zmax )
            error[elno] = maximum(abs.( Array(intrp_brck.v[st:en])-fex[:]))
println("$elno). error = ", error[elno])
        else
println("$elno). No interpolation points in this element")
        end
    end

#    println("==============================================")
#    println("l_infinity interpolation error in each element")
#    display(error)
    l_infinity_local = maximum(error)
    l_infinity_domain = MPI.Allreduce(l_infinity_local, MPI.MAX, mpicomm)
    println("First run_brick_interpolation_test(): l_infinity interpolation error in domain")
    display(l_infinity_domain)
#    pid = MPI.Comm_rank(mpicomm)
#    npr = MPI.Comm_size(mpicomm)

#    for i in 0:npr-1
#        if i == pid
#            println("pid = $pid; l_infinity_local = $l_infinity_local; l_infinity_domain = $l_infinity_domain")
#        end
#        MPI.Barrier(mpicomm)
#    end

#    return l_infinity_domain < 1.0e-14
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
#Base.@kwdef struct TestSphereSetup{FT}
#  p_ground::FT = MSLP
#  T_initial::FT = 255
#  domain_height::FT = 30e3
#end


struct TestSphereSetup{DT}
  p_ground::DT 
  T_initial::DT 
  domain_height::DT
  
  function TestSphereSetup(p_ground::DT, T_initial::DT, domain_height::DT) where DT <: AbstractFloat
    return new{DT}(p_ground, T_initial, domain_height)
  end
end
#----------------------------------------------------------------------------
# Cubed sphere, lat/long interpolation test
#----------------------------------------------------------------------------
function run_cubed_sphere_interpolation_test()
    CLIMA.init()

    FT = Float64 #Float32 #Float64
    mpicomm = MPI.COMM_WORLD
    root = 0

    ll = uppercase(get(ENV, "JULIA_LOG_LEVEL", "INFO"))
    loglevel = Dict("DEBUG" => Logging.Debug,
                    "WARN"  => Logging.Warn,
                    "ERROR" => Logging.Error,
                    "INFO"  => Logging.Info)[ll]

    logger_stream = MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull
    global_logger(ConsoleLogger(logger_stream, loglevel))

    domain_height = FT(30e3) 

    polynomialorder = 5#12#5# 12#1#4 #5
    numelem_horz = 3#4 #6
    numelem_vert = 4#1 #1 #1#6 #8

    #-------------------------
    _x, _y, _z = CLIMA.Mesh.Grids.vgeoid.x1id, CLIMA.Mesh.Grids.vgeoid.x2id, CLIMA.Mesh.Grids.vgeoid.x3id
    _ρ, _ρu, _ρv, _ρw = 1, 2, 3, 4
    #-------------------------
    vert_range = grid1d(FT(planet_radius), FT(planet_radius + domain_height), nelem = numelem_vert)

 #   vert_range = grid1d(FT(1.0), FT(2.0), nelem = numelem_vert)
  
    lat_res  = FT( 1 * π / 180.0) # 5 degree resolution
    long_res = FT( 1 * π / 180.0) # 5 degree resolution
    nel_grd  = 10
    r_res    = FT((vert_range[end] - vert_range[1])/FT(nel_grd)) #1000.00    # 1000 m vertical resolution

    #----------------------------------------------------------
    setup = TestSphereSetup(FT(MSLP),FT(255),FT(30e3))

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
                       NoSubsidence{FT}(),
                       nothing, 
                       NoFluxBC(),
                       setup)

    dg = DGModel(model, grid, Rusanov(),
                 CentralNumericalFluxDiffusive(), CentralGradPenalty())

    Q = init_ode_state(dg, FT(0))

    device = typeof(Q.data) <: Array ? CPU() : CUDA()
    #------------------------------
    x1 = @view grid.vgeo[:,_x,:]
    x2 = @view grid.vgeo[:,_y,:]
    x3 = @view grid.vgeo[:,_z,:]

    xmax = maximum( abs.(x1) )
    ymax = maximum( abs.(x2) )
    zmax = maximum( abs.(x3) )

    st_idx = _ρ # state vector

    var = @view Q.data[:,st_idx,:]

#    fcn(x,y,z) = x .* y .* z # sample function
    fcn(x,y,z) = sin.(x) .* cos.(y) .* cos.(z) # sample function

    var .= fcn( x1 ./ xmax, x2 ./ ymax, x3 ./ zmax )
  #------------------------------
    qm1 = polynomialorder + 1

    intrp_cs = InterpolationCubedSphere(grid, collect(vert_range), numelem_horz, lat_res, long_res, r_res)
println("After setting up intrp_cs")
    interpolate_cubed_sphere!(intrp_cs.offset, intrp_cs.m1_r, intrp_cs.wb, intrp_cs.ξ1, intrp_cs.ξ2, intrp_cs.ξ3, intrp_cs.flg, intrp_cs.fac, intrp_cs.v, Q.data, st_idx)#,Val(qm1))
    for i in 1:20
        @time interpolate_cubed_sphere!(intrp_cs.offset, intrp_cs.m1_r, intrp_cs.wb, intrp_cs.ξ1, intrp_cs.ξ2, intrp_cs.ξ3, intrp_cs.flg, intrp_cs.fac, intrp_cs.v, Q.data, st_idx)#,Val(qm1))
    end
#    @launch(device, threads=(qm1,qm1), blocks=1, interpolate_cubed_sphere!(intrp_cs, Q.data, st_idx, Val(qm1)))
    #----------------------------------------------------------
    Nel = length( grid.topology.realelems )

    error = zeros(FT, Nel) 
        
    offset = Array(intrp_cs.offset)
    radc   = Array(intrp_cs.radc)
    latc   = Array(intrp_cs.latc) 
    longc  = Array(intrp_cs.longc) 
    v      = Array(intrp_cs.v)

    for elno in 1:Nel
        st  = offset[elno] + 1
        en  = offset[elno+1]
        if en ≥ st 

            x1_grd = radc[st:en] .* sin.(latc[st:en]) .* cos.(longc[st:en]) # inclination -> latitude; azimuthal -> longitude.
            x2_grd = radc[st:en] .* sin.(latc[st:en]) .* sin.(longc[st:en]) # inclination -> latitude; azimuthal -> longitude.
            x3_grd = radc[st:en] .* cos.(latc[st:en])
        
            fex = fcn( x1_grd ./ xmax , x2_grd ./ ymax , x3_grd ./ zmax )
            error[elno] = maximum(abs.(v[st:en] - fex[:]))
        end
    end
    #----------------------------------------------------------
    l_infinity_local = maximum(error)
    l_infinity_domain = MPI.Allreduce(l_infinity_local, MPI.MAX, mpicomm)
    toler = 2.0e-5
    println("l_infinity_domain = $l_infinity_domain; toler = $toler")

#    if l_infinity_domain > toler
#        pid = MPI.Comm_rank(mpicomm)
#        npr = MPI.Comm_size(mpicomm)

#        for i in 0:npr-1
#            if i == pid
#                println("pid = $pid; l_infinity_local = $l_infinity_local; l_infinity_domain = $l_infinity_domain")
#            end
#            MPI.Barrier(mpicomm)
#        end
        
#        if pid == 0
#            println("l_infinity_domain = $l_infinity_domain; toler = $toler")
#        end
#    end
#----------------------------------------------------------------------------
    return l_infinity_domain < toler # 1.0e-12
return true
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

#@testset "Interpolation tests" begin
#    @test run_brick_interpolation_test()
#    @test run_cubed_sphere_interpolation_test()
#end
#run_brick_interpolation_test()
run_cubed_sphere_interpolation_test()
#------------------------------------------------

