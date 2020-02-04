using Test, MPI
import GaussQuadrature
using CLIMA
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Mesh.Geometry
using CLIMA.Mesh.Interpolation
using StaticArrays
using GPUifyLoops

using CLIMA.VariableTemplates
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
#------------------------------------------------
#if !@isdefined integration_testing
#    const integration_testing =
#    parse(Bool, lowercase(get(ENV,"JULIA_CLIMA_INTEGRATION_TESTING","false")))
#end
#------------------------------------------------
function run_brick_interpolation_test()
  MPI.Initialized() || MPI.Init()
  DA = CLIMA.array_type()
#@testset "LocalGeometry" begin
    FT = Float64
#    ArrayType = ArrayType #Array
    mpicomm = MPI.COMM_WORLD
    pid = MPI.Comm_rank(mpicomm)
    npr = MPI.Comm_size(mpicomm)

    xmin, ymin, zmin = 0, 0, 0                   # defining domain extent
    xmax, ymax, zmax = 2000, 400, 2000
#    xres = [FT(200), FT(200), FT(200)] # resolution of interpolation grid
#    xres = [FT(400), FT(200), FT(400)] # resolution of interpolation grid
#    xres = [FT(5), FT(5), FT(5)] # resolution of interpolation grid
    xres = [FT(10), FT(10), FT(10)] # resolution of interpolation grid

#    Ne        = (4,4,4)
    Ne        = (20,4,20)
#    Ne        = (40,8,40)
#    Ne        = (40,16,40)
#    Ne        = (60,32,60)
#    Ne        = (80,64,80)
#    Ne        = (4,2,4)

    polynomialorder = 5 #8# 5 #8 #4
    #-------------------------
    _x, _y, _z = CLIMA.Mesh.Grids.vgeoid.x1id, CLIMA.Mesh.Grids.vgeoid.x2id, CLIMA.Mesh.Grids.vgeoid.x3id
    _ρ, _ρu, _ρv, _ρw = 1, 2, 3, 4
    #-------------------------

    brickrange = (range(FT(xmin); length=Ne[1]+1, stop=xmax),
                  range(FT(ymin); length=Ne[2]+1, stop=ymax),
                  range(FT(zmin); length=Ne[3]+1, stop=zmax))
    topl = StackedBrickTopology(mpicomm, brickrange, periodicity = (true, true, false))
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
 

    #fcn(x,y,z) = x .* y .* z # sample function
    fcn(x,y,z) = sin.(x) .* cos.(y) .* cos.(z) # sample function

    

    #----calling interpolation function on state variable # st_idx--------------------------
    nvars = size(Q.data,2)

    for vari in 1:nvars
		Q.data[:,vari,:] = fcn( x1 ./ xmax, x2 ./ ymax, x3 ./ zmax )
    end

    xbnd = Array{FT}(undef,2,3)
 
    xbnd[1,1] = FT(xmin); xbnd[2,1] = FT(xmax)
    xbnd[1,2] = FT(ymin); xbnd[2,2] = FT(ymax)
    xbnd[1,3] = FT(zmin); xbnd[2,3] = FT(zmax)

    if pid==0 
        println("interpolation setup")
    end
    intrp_brck = InterpolationBrick(grid, xbnd, xres)
#    np_ig = Array(intrp_brck.offset)[end]    
#    iv = DA( Array{FT}(undef,np_ig,nvars) )
 
    if pid==0 
    	println("timing interpolation/ first call")
    end
    svi = interpolate_brick!(intrp_brck, Q.data)

#    println("timing interpolation/ multiple calls")
#    for i in 1:10
#    	@time svi = interpolate_brick!(intrp_brck, Q.data)
#    end

    err_inf_dom = zeros(FT, nvars) 
        
    x1g = intrp_brck.x1g
    x2g = intrp_brck.x2g
    x3g = intrp_brck.x3g


    if pid==0
        nx1 = length(x1g); nx2 = length(x2g); nx3 = length(x3g)
		x1 = Array{FT}(undef,nx1,nx2,nx3); x2 = similar(x1); x3 = similar(x1)

        for k in 1:nx3, j in 1:nx2, i in 1:nx1
            x1[i,j,k] = x1g[i] 
   	        x2[i,j,k] = x2g[j]  
       	    x3[i,j,k] = x3g[k]
        end
		fex = fcn( x1 ./ xmax , x2 ./ ymax , x3 ./ zmax )

        for vari in 1:nvars
			err_inf_dom[vari] = maximum( abs.(svi[:,:,:,vari] .- fex[:,:,:]) )
        end
    end
#    println("timing interpolate_brick")
#    for i in 1:10
#        @time interpolate_brick_local!(intrp_brck, Q.data, iv)
#    end

	toler = 1.0E-9
    return maximum(err_inf_dom) < toler #1.0e-14
#return true
    #----------------
end #function run_brick_interpolation_test


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
# thermodynamic variables of interest
function test_vars_thermo(FT)
  @vars begin
    q_liq::FT
    q_ice::FT
    q_vap::FT
    T::FT
    θ_liq_ice::FT
    θ_dry::FT
    θ_v::FT
    e_int::FT
    h_m::FT
    h_t::FT
  end
end
#----------------------------------------------------------------------------
# Cubed sphere, lat/long interpolation test
#----------------------------------------------------------------------------
function run_cubed_sphere_interpolation_test()
    CLIMA.init()

    DA = CLIMA.array_type()
    FT = Float64 #Float32 #Float64
    mpicomm = MPI.COMM_WORLD
    root = 0
    pid  = MPI.Comm_rank(mpicomm)
    npr  = MPI.Comm_size(mpicomm)


    ll = uppercase(get(ENV, "JULIA_LOG_LEVEL", "INFO"))
    loglevel = Dict("DEBUG" => Logging.Debug,
                    "WARN"  => Logging.Warn,
                    "ERROR" => Logging.Error,
                    "INFO"  => Logging.Info)[ll]

    logger_stream = MPI.Comm_rank(mpicomm) == 0 ? stderr : devnull
    global_logger(ConsoleLogger(logger_stream, loglevel))

    domain_height = FT(30e3) 

    polynomialorder = 5#5#12#5# 12#1#4 #5
    numelem_horz = 6#4 #6
    numelem_vert = 4 #4#1 #1 #1#6 #8

    #-------------------------
    _x, _y, _z = CLIMA.Mesh.Grids.vgeoid.x1id, CLIMA.Mesh.Grids.vgeoid.x2id, CLIMA.Mesh.Grids.vgeoid.x3id
    _ρ, _ρu, _ρv, _ρw = 1, 2, 3, 4
    #-------------------------
    vert_range = grid1d(FT(planet_radius), FT(planet_radius + domain_height), nelem = numelem_vert)

 #   vert_range = grid1d(FT(1.0), FT(2.0), nelem = numelem_vert)
  
#    lat_res  = FT( 10 * π / 180.0) # 5 degree resolution
#    long_res = FT( 10 * π / 180.0) # 5 degree resolution
    #lat_res  = FT( 1 * π / 180.0) # 5 degree resolution
    #long_res = FT( 1 * π / 180.0) # 5 degree resolution
    lat_res  = FT( 1.0 * π / 180.0) # 5 degree resolution
    long_res = FT( 1.0 * π / 180.0) # 5 degree resolution
    #nel_vert_grd  = 100 #100 #50 #10#50
    nel_vert_grd  = 20 #100 #50 #10#50
    r_res    = FT((vert_range[end] - vert_range[1])/FT(nel_vert_grd)) #1000.00    # 1000 m vertical resolution

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

    xmax = FT(planet_radius) #maximum( abs.(x1) )
    ymax = FT(planet_radius) #maximum( abs.(x2) )
    zmax = FT(planet_radius) #maximum( abs.(x3) )

    #fcn(x,y,z) = sin.(x) .* cos.(y) .* cos.(z) # sample function
    fcn(x,y,z) = x .* y .* z # sample function

    nvars = size(Q.data,2)

    for i in 1:nvars
		Q.data[:,i,:] .= fcn( x1 ./ xmax, x2 ./ ymax, x3 ./ zmax )
    end
  #------------------------------
    qm1 = polynomialorder + 1

    intrp_cs = InterpolationCubedSphere(grid, collect(vert_range), numelem_horz, lat_res, long_res, r_res);
#    np_ig = tmp[end]
#for i in 1:10
#    interp_cs = []

#    @time intrp_cs = InterpolationCubedSphere(grid, collect(vert_range), numelem_horz, lat_res, long_res, r_res)
#end

    svi = interpolate_cubed_sphere!(intrp_cs, Q.data) # interpolates state vector onto svi
#    for i in 1:20
#        @time interpolate_cubed_sphere_local!(intrp_cs, Q.data, iv)#,Val(qm1))
#    end
    #----------------------------------------------------------
#    Nel = length( grid.topology.realelems )

    err_inf_dom = zeros(FT, nvars) 
        
    rad   = intrp_cs.rad_grd
    lat   = intrp_cs.lat_grd
    long  = intrp_cs.long_grd

    if pid==0
        nrad = length(rad); nlat = length(lat); nlong = length(long)
		x1g = Array{FT}(undef,nrad,nlat,nlong); x2g = similar(x1g); x3g = similar(x1g)

        for k in 1:nlong, j in 1:nlat, i in 1:nrad
            x1g[i,j,k] = rad[i] * sin(lat[j]) * cos(long[k]) # inclination -> latitude; azimuthal -> longitude.
   	        x2g[i,j,k] = rad[i] * sin(lat[j]) * sin(long[k]) # inclination -> latitude; azimuthal -> longitude.
       	    x3g[i,j,k] = rad[i] * cos(lat[j])
        end
		fex = fcn( x1g ./ xmax , x2g ./ ymax , x3g ./ zmax )

        for vari in 1:nvars
			err_inf_dom[vari] = maximum( abs.(svi[:,:,:,vari] .- fex[:,:,:]) )
        end
    end

	MPI.Bcast!(err_inf_dom, root, mpicomm)

    toler = 1.0E-7
    return maximum(err_inf_dom) < toler # 1.0e-12
end 
#----------------------------------------------------------------------------

@testset "Interpolation tests" begin
    @test run_brick_interpolation_test()
    @test run_cubed_sphere_interpolation_test()
end
#run_brick_interpolation_test()
#run_cubed_sphere_interpolation_test()
#------------------------------------------------

