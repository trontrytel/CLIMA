using CLIMA.PlanetParameters
export PeriodicBC, NoFluxBC, InitStateBC, DYCOMS_BC, RayleighBenardBC, SurfaceDrivenBubbleBC

function atmos_boundary_flux_diffusive!(nf::CentralNumericalFluxDiffusive, bc,
                                        atmos::AtmosModel, 
                                        F⁺,Y⁺, Σ⁺,α⁺, 
                                        n⁻, 
                                        F⁻,Y⁻,Σ⁻,α⁻,
                                        bctype, t, 
                                        Y₁⁻, Σ₁⁻, α₁⁻)
  # Y  ≡ state vars
  # α  ≡ auxiliary vars
  # Σ  ≡ diffusive vars
  # F  ≡ flux
  # X₁ ≡ X at the first interior node
  # bctype ≡ `wall` identifier
  # t ≡ simulation time
  
  FT = eltype(F⁺)
  atmos_boundary_state!(nf, bc, atmos, Y⁺, Σ⁺, α⁺, n⁻,
                        Y⁻, Σ⁻, α⁻, bctype, t,
                        Y₁⁻, Σ₁⁻, α₁⁻)
  fill!(parent(F⁺), -zero(FT))
  flux_diffusive!(atmos, F⁺, Y⁺, Σ⁺, α⁺, t)
end

#TODO: figure out a better interface for this.
# at the moment we can just pass a function, but we should do something better
# need to figure out how subcomponents will interact.
function atmos_boundary_state!(::Rusanov, f::Function, m::AtmosModel,
                               Y⁺::Vars, α⁺::Vars, 
                               n⁻, 
                               Y⁻::Vars, α⁻::Vars, 
                               bctype, t, _...)
  f(Y⁺, α⁺, n⁻, Y⁻, α⁻, bctype, t)
end

function atmos_boundary_state!(::CentralNumericalFluxDiffusive, f::Function,
                               m::AtmosModel, 
                               Y⁺::Vars, Σ⁺::Vars, α⁺::Vars, 
                               n⁻, 
                               Y⁻::Vars, Σ⁻::Vars,α⁻::Vars, 
                               bctype, t, _...)
  f(Y⁺, Σ⁺, α⁺, n⁻, Y⁻, Σ⁻, α⁻, bctype, t)
end

# lookup boundary condition by face
function atmos_boundary_state!(nf::Rusanov, bctup::Tuple, m::AtmosModel,
                               Y⁺::Vars, α⁺::Vars, n⁻, Y⁻::Vars,
                               α⁻::Vars, bctype, t, _...)
  atmos_boundary_state!(nf, bctup[bctype], m, Y⁺, α⁺, n⁻, Y⁻, α⁻,
                        bctype, t)
end

function atmos_boundary_state!(nf::CentralNumericalFluxDiffusive,
                               bctup::Tuple, m::AtmosModel, Y⁺::Vars,
                               Σ⁺::Vars, α⁺::Vars, n⁻, Y⁻::Vars,
                               Σ⁻::Vars, α⁻::Vars, bctype, t, _...)
  atmos_boundary_state!(nf, bctup[bctype], m, Y⁺, Σ⁺, α⁺, n⁻, Y⁻,
                        Σ⁻, α⁻, bctype, t)
end


abstract type BoundaryCondition
end

"""
    PeriodicBC <: BoundaryCondition

Assume that the topology is periodic and hence nothing special needs to be done at the boundaries.
"""
struct PeriodicBC <: BoundaryCondition end

# TODO: assert somewhere that the topology is actually periodic when using those
atmos_boundary_state!(_, ::PeriodicBC, _...) = nothing

"""
    NoFluxBC <: BoundaryCondition

Set the momentum at the boundary to be zero.

# TODO: This should be fixed later once BCs are figured out (likely want
# different things here?)

"""
struct NoFluxBC <: BoundaryCondition
end

function atmos_boundary_state!(::Rusanov, bc::NoFluxBC, m::AtmosModel,
                               Y⁺::Vars, α⁺::Vars, n⁻, Y⁻::Vars,
                               α⁻::Vars, bctype, t, _...)
  FT = eltype(Y⁻)
  Y⁺.ρ = Y⁻.ρ
  Y⁺.ρu -= 2 * dot(Y⁻.ρu, n⁻) * SVector(n⁻)
end

function atmos_boundary_state!(::CentralNumericalFluxDiffusive, bc::NoFluxBC,
                               m::AtmosModel, Y⁺::Vars, Σ⁺::Vars,
                               α⁺::Vars, n⁻, Y⁻::Vars, Σ⁻::Vars,
                               α⁻::Vars, bctype, t, _...)
  FT = eltype(Y⁻)
  Y⁺.ρ = Y⁻.ρ
  Y⁺.ρu -= 2 * dot(Y⁻.ρu, n⁻) * SVector(n⁻)
  
  fill!(getfield(Σ⁺, :array), FT(0))
end

"""
    InitStateBC <: BoundaryCondition

Set the value at the boundary to match the `init_state!` function. This is
mainly useful for cases where the problem has an explicit solution.

# TODO: This should be fixed later once BCs are figured out (likely want
# different things here?)
"""
struct InitStateBC <: BoundaryCondition
end
function atmos_boundary_state!(::Rusanov, bc::InitStateBC, m::AtmosModel,
                               Y⁺::Vars, α⁺::Vars, n⁻, Y⁻::Vars,
                               α⁻::Vars, bctype, t, _...)
  init_state!(m, Y⁺, α⁺, α⁺.coord, t)
end
function atmos_boundary_state!(::CentralNumericalFluxDiffusive, bc::InitStateBC,
                               m::AtmosModel, Y⁺::Vars, Σ⁺::Vars,
                               α⁺::Vars, n⁻, Y⁻::Vars, Σ⁻::Vars,
                               α⁻::Vars, bctype, t, _...)
  init_state!(m, Y⁺, α⁺, α⁺.coord, t)
end


"""
  DYCOMS_BC <: BoundaryCondition
  Prescribes boundary conditions for Dynamics of Marine Stratocumulus Case
"""
struct DYCOMS_BC{FT} <: BoundaryCondition
  C_drag::FT
  LHF::FT
  SHF::FT
end
function atmos_boundary_state!(::Rusanov, bc::DYCOMS_BC, m::AtmosModel,
                               Y⁺::Vars, α⁺::Vars, n⁻, Y⁻::Vars,
                               α⁻::Vars, bctype, t, Y₁⁻::Vars, α₁::Vars)
  # Y⁻ is the 𝐘⁻ state while Y⁺ is the 𝐘⁺ state at an interface.
  # at the boundaries the ⁻, minus side states are the interior values
  # Y₁⁻ is 𝐘 at the first interior nodes relative to the bottom wall
  FT = eltype(Y⁺)
  # Get values from minus-side state
  ρ⁻ = Y⁻.ρ
  UM, VM, WM = Y⁻.ρu
  EM = Y⁻.ρe
  QTM = Y⁻.moisture.ρq_tot
  uM, vM, wM  = UM/ρ⁻, VM/ρ⁻, WM/ρ⁻
  q_totM = QTM/ρ⁻
  Un⁻ = n⁻[1] * UM + n⁻[2] * VM + n⁻[3] * WM

  # Assign reflection wall boundaries (top wall)
  Y⁺.ρu = SVector(UM - 2 * n⁻[1] * Un⁻,
                      VM - 2 * n⁻[2] * Un⁻,
                      WM - 2 * n⁻[3] * Un⁻)

  # Assign scalar values at the boundaries
  Y⁺.ρ = ρ⁻
  Y⁺.moisture.ρq_tot = QTM
end
function atmos_boundary_flux_diffusive!(nf::CentralNumericalFluxDiffusive,
                                        bc::DYCOMS_BC, atmos::AtmosModel,
                                        F⁺, Y⁺, Σ⁺, α⁺, n⁻,
                                        F⁻, Y⁻, Σ⁻, α⁻,
                                        bctype, t,
                                        Y₁⁻, Σ₁⁻, α₁⁻)
  FT = eltype(Y⁺)

  # Y⁻ is the 𝐘⁻ state while Y⁺ is the 𝐘⁺ state at an interface.
  # at the boundaries the ⁻, minus side states are the interior values
  # Y₁⁻⁻ is 𝐘 at the first interior nodes relative to the bottom wall
  # Get values from minus-side state
  ρ⁻ = Y⁻.ρ
  U⁻, V⁻, W⁻ = Y⁻.ρu
  E⁻ = Y⁻.ρe
  QT⁻ = Y⁻.moisture.ρq_tot
  u⁻, v⁻, w⁻  = U⁻/ρ⁻, V⁻/ρ⁻, W⁻/ρ⁻
  q_tot⁻ = QT⁻/ρ⁻
  Un⁻ = n⁻[1] * U⁻ + n⁻[2] * V⁻ + n⁻[3] * W⁻

  # Assign reflection wall boundaries (top wall)
  Y⁺.ρu = SVector(U⁻ - 2 * n⁻[1] * Un⁻,
                      V⁻ - 2 * n⁻[2] * Un⁻,
                      W⁻ - 2 * n⁻[3] * Un⁻)

  # Assign scalar values at the boundaries
  Y⁺.ρ = ρ⁻
  Y⁺.moisture.ρq_tot = QT⁻
  # Assign diffusive fluxes at boundaries
  Σ⁺ = Σ⁻
  if bctype != 1
    flux_diffusive!(atmos, F⁺, Y⁺, Σ⁺, α⁺, t)
  else
    # ------------------------------------------------------------------------
    # (<var>_FN) First node values (First interior node from bottom wall)
    # ------------------------------------------------------------------------
    z_FN             = α₁⁻.coord[3]
    ρ_FN             = Y₁⁻.ρ
    U_FN, V_FN, W_FN = Y₁⁻.ρu
    E_FN             = Y₁⁻.ρe
    u_FN, v_FN, w_FN = U_FN/ρ_FN, V_FN/ρ_FN, W_FN/ρ_FN
    windspeed_FN     = sqrt(u_FN^2 + v_FN^2 + w_FN^2)
    q_tot_FN         = Y₁⁻.moisture.ρq_tot / ρ_FN
    e_int_FN         = E_FN/ρ_FN - windspeed_FN^2/2 - grav*z_FN
    TS_FN            = PhaseEquil(e_int_FN, ρ_FN, q_tot_FN)
    T_FN             = air_temperature(TS_FN)
    q_vap_FN         = q_tot_FN - PhasePartition(TS_FN).liq
    # --------------------------
    # Bottom boundary quantities
    # --------------------------
    z⁻          = α⁻.coord[3]
    q_tot⁻      = QT⁻/ρ⁻
    windspeed   = sqrt(u⁻^2 + v⁻^2 + w⁻^2)
    e_int⁻      = E⁻/ρ⁻ - windspeed^2/2 - grav*z⁻
    TS⁻         = PhaseEquil(e_int⁻, ρ⁻, q_tot⁻)
    q_vap⁻      = q_tot⁻ - PhasePartition(TS⁻).liq
    T⁻          = air_temperature(TS⁻)
    # ----------------------------------------------------------
    # Extract components of diffusive momentum flux (minus-side)
    # ----------------------------------------------------------
    _, τ⁻ = turbulence_tensors(atmos.turbulence, Y⁻, Σ⁻, α⁻, t)

    # ----------------------------------------------------------
    # Boundary momentum fluxes
    # ----------------------------------------------------------
    # Case specific for flat bottom topography, normal vector is n⃗ = k⃗ = [0, 0, 1]ᵀ
    # A more general implementation requires (n⃗ ⋅ ∇A) to be defined where A is replaced by the appropriate flux terms
    C_drag = bc.C_drag
    τ13⁺  = - C_drag * windspeed_FN * u_FN
    τ23⁺  = - C_drag * windspeed_FN * v_FN
    # Assign diffusive momentum and moisture fluxes
    # (i.e. ρ𝛕 terms)
    τ⁺ = SHermitianCompact{3, FT, 6}(SVector(FT(0), τ⁻[2,1], τ13⁺, FT(0), τ23⁺,
                                             FT(0)))

    # ----------------------------------------------------------
    # Boundary moisture fluxes
    # ----------------------------------------------------------
    # really ∇q_tot is being used to store d_q_tot
    d_q_tot⁺  = SVector(FT(0), FT(0), bc.LHF/(LH_v0))

    # ----------------------------------------------------------
    # Boundary energy fluxes
    # ----------------------------------------------------------
    # Assign diffusive enthalpy flux (i.e. ρ(J+D) terms)
    d_h_tot⁺ = SVector(FT(0), FT(0), bc.LHF + bc.SHF)

    flux_diffusive!(atmos, F⁺, Y⁺, τ⁺, d_h_tot⁺)
    flux_diffusive!(atmos.moisture, F⁺, Y⁺, d_q_tot⁺)
  end
end

"""
  RayleighBenardBC <: BoundaryCondition

# Fields
$(DocStringExtensions.FIELDS)
"""
struct RayleighBenardBC{FT} <: BoundaryCondition
  "Prescribed bottom wall temperature [K]"
  T_bot::FT
  "Prescribed top wall temperature [K]"
  T_top::FT
end
# Rayleigh-Benard problem with two fixed walls (prescribed temperatures)
function atmos_boundary_state!(::Rusanov, bc::RayleighBenardBC, m::AtmosModel,
                               Y⁺::Vars, α⁺::Vars, n⁻, Y⁻::Vars,
                               α⁻::Vars, bctype, t,_...)
  # Dry Rayleigh Benard Convection
  @inbounds begin
    FT = eltype(Y⁺)
    Y⁺.ρu = 2 * dot(Y⁻.ρu, n⁻) * SVector(n⁻) - Y⁻.ρu
    if bctype == 1
      ρe_int⁺ = Y⁺.ρ * cv_d * (bc.T_bot - T_0)
    else
      ρe_int⁺ = Y⁺.ρ * cv_d * (bc.T_top - T_0)
    end
    Y⁺.ρe = 2 * (ρe_int⁺ + Y⁺.ρ * α⁺.coord[3] * grav) - Y⁻.ρe
    nothing
  end
end
function atmos_boundary_state!(::CentralNumericalFluxDiffusive, bc::RayleighBenardBC,
                               m::AtmosModel, Y⁺::Vars, Σ⁺::Vars,
                               α⁺::Vars, n⁻, Y⁻::Vars, Σ⁻::Vars,
                               α⁻::Vars, bctype, t, _...)
  # Dry Rayleigh Benard Convection
  @inbounds begin
    FT = eltype(Y⁻)
    Y⁺.ρu = 2 * dot(Y⁻.ρu, n⁻) * SVector(n⁻) - Y⁻.ρu
    if bctype == 1
      ρe_int⁺ = Y⁺.ρ * cv_d * (bc.T_bot - T_0)
    else
      ρe_int⁺ = Y⁺.ρ * cv_d * (bc.T_top - T_0)
    end
    Y⁺.ρe = (ρe_int⁺ + Y⁺.ρ * α⁺.coord[3] * grav)
    nothing
  end
end


"""
  SurfaceDrivenBubbleBC <: BoundaryCondition

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SurfaceDrivenBubbleBC{FT} <: BoundaryCondition
  "Prescribed MSEF Magnitude [W/m^2]"
  F₀::FT
  "Spatial Parameter [m]"
  σ::FT
  "Surface Heater Radius [m]"
  a::FT
  "Surface Heater Center [m]"
  x₀::FT
  "Time Cutoff [s]"
  t₁::FT
end
function atmos_boundary_state!(::Rusanov, bc::SurfaceDrivenBubbleBC, 
                               m::AtmosModel,
                               Y⁺::Vars, 
                               α⁺::Vars, 
                               n⁻, 
                               Y⁻::Vars,
                               α⁻::Vars, 
                               bctype, t,_...)
  FT = eltype(Y⁺)
  @inbounds begin
    # Momentum b.c. prescribed (no flow across wall boundary)
    Y_bc  = dot(Y⁻.ρu, n⁻)*SVector(n⁻)
    Y⁺.ρu = -Y⁻.ρu + 2*Y_bc
  end
  nothing
end
function atmos_boundary_state!(::CentralNumericalFluxDiffusive, bc::SurfaceDrivenBubbleBC,
                               m::AtmosModel, 
                               Y⁺::Vars, 
                               Σ⁺::Vars,
                               α⁺::Vars, 
                               n⁻, 
                               Y⁻::Vars, 
                               Σ⁻::Vars,
                               α⁻::Vars, 
                               bctype, t, _...)
  FT = eltype(Y⁻)
  k̂  = α⁻.orientation.∇Φ / norm(α⁻.orientation.∇Φ)  
  r = sqrt((α⁻.coord[1]-bc.x₀)^2 + (α⁻.coord[2]-bc.x₀)^2) 
  F₀ =  bc.F₀*(1 - sign(t-bc.t₁))/2 
  @inbounds begin
    # Energy flux prescribed (diffusive flux through bottom wall)
    # MSEF ≡ Moist Static Energy Flux
    if bctype == 1
      if r > bc.a 
        MSEF      = F₀ * exp(-(r-bc.a)^2 / bc.σ^2)
        Σ⁺.∇h_tot = -Σ⁻.∇h_tot + 2 * MSEF * k̂
      else
        MSEF      = F₀
        Σ⁺.∇h_tot = -Σ⁻.∇h_tot + 2 * MSEF * k̂
      end
    end
  end
  nothing
end

