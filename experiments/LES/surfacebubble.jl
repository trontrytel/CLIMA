using Distributions
using Random
using StaticArrays
using Test

using CLIMA
using CLIMA.Atmos
using CLIMA.GenericCallbacks
using CLIMA.LowStorageRungeKuttaMethod
using CLIMA.Mesh.Filters
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters
using CLIMA.VariableTemplates


"""
  Surface Driven Thermal Bubble
"""
function init_surfacebubble!(state, aux, (x,y,z), t)
  FT            = eltype(state)
  R_gas::FT     = R_d
  c_p::FT       = cp_d
  c_v::FT       = cv_d
  γ::FT         = c_p / c_v
  p0::FT        = MSLP

  xc::FT        = 1250
  yc::FT        = 1250
  zc::FT        = 1000
  r             = sqrt((x-xc)^2+(y-yc)^2+(z-zc)^2)
  rc::FT        = 500
  θ_ref::FT     = 300
  Δθ::FT        = 0

  if r <= rc
    Δθ          = FT(5) * cospi(r/rc/2)
  end

  #Perturbed state:
  θ            = θ_ref + Δθ # potential temperature
  π_exner      = FT(1) - grav / (c_p * θ) * z # exner pressure
  ρ            = p0 / (R_gas * θ) * (π_exner)^ (c_v / R_gas) # density
  P            = p0 * (R_gas * (ρ * θ) / p0) ^(c_p/c_v) # pressure (absolute)
  T            = P / (ρ * R_gas) # temperature
  ρu           = SVector(FT(0),FT(0),FT(0))
  # energy definitions
  e_kin        = FT(0)
  e_pot        = grav * z
  ρe_tot       = ρ * total_energy(e_kin, e_pot, T)
  state.ρ      = ρ
  state.ρu     = ρu
  state.ρe     = ρe_tot
  state.moisture.ρq_tot = FT(0)
end

function config_surfacebubble(FT, N, resolution, xmax, ymax, zmax)
    
  #
    # Sponge
    # 
    c_sponge = 1
    zsponge = FT(1500.0)
    u_relaxation = SVector{3,FT}(0,0,0)
    rayleigh_sponge = RayleighSponge{FT}(zmax, zsponge, c_sponge, u_relaxation, 2)
    
    #
    # Boundary conditions
    #
    F₀ = FT(150)
    σ = FT(100)
    a = FT(50)
    t₁ = FT(500)
    x₀ = FT(500)
    bc = SurfaceDrivenBubbleBC{FT}(F₀, σ, a, x₀, t₁)

    C_smag = FT(0.23)
    
    implicitsolver = CLIMA.DefaultSolverType()

    config = CLIMA.LES_Configuration("SurfaceDrivenBubble", 
                                     N, resolution, xmax, ymax, zmax,
                                     init_surfacebubble!,
                                     solver_type=implicitsolver,
                                     C_smag=C_smag,
                                     moisture=EquilMoist(),
                                     sources=(Gravity(), rayleigh_sponge),
                                     bc=bc)
    return config
end

function main()
    CLIMA.init()
    FT = Float64
    # DG polynomial order
    N = 4
    # Domain resolution and size
    Δh = FT(50)
    Δv = FT(50)
    resolution = (Δh, Δh, Δv)
    xmax = 2500
    ymax = 2500
    zmax = 5000
    t0 = FT(0)
    timeend = FT(2000)
    
    CFL_max = FT(0.1)

    driver_config = config_surfacebubble(FT, N, resolution, xmax, ymax, zmax)
    solver_config = CLIMA.setup_solver(t0, timeend, Courant_number=CFL_max, driver_config, forcecpu=true)

    cbtmarfilter = GenericCallbacks.EveryXSimulationSteps(1) do (init=false)
        Filters.apply!(solver_config.Q, 6, solver_config.dg.grid, TMARFilter())
        nothing
    end

    result = CLIMA.invoke!(solver_config;
                          user_callbacks=(cbtmarfilter,),
                          check_euclidean_distance=true)
end

main()
