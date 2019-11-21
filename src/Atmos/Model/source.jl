using CLIMA.PlanetParameters: Omega
export Gravity, RayleighSponge, Subsidence, GeostrophicForcing, Coriolis, BulkDrag

# kept for compatibility
# can be removed if no functions are using this
function atmos_source!(f::Function, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  f(source, state, aux, t)
end
function atmos_source!(::Nothing, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
end
# sources are applied additively
function atmos_source!(stuple::Tuple, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  map(s -> atmos_source!(s, m, source, state, aux, t), stuple)
end

abstract type Source
end

struct Gravity <: Source
end
function atmos_source!(::Gravity, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  source.ρu -= state.ρ * aux.orientation.∇Φ
end

struct Subsidence <: Source
end
function atmos_source!(::Subsidence, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  n = aux.orientation.∇Φ ./ norm(aux.orientation.∇Φ)
  source.ρu -= m.radiation.D_subsidence * dot(state.ρu, n) * n
end

struct Coriolis <: Source
end
function atmos_source!(::Coriolis, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  # note: this assumes a SphericalOrientation
  source.ρu -= SVector(0, 0, 2*Omega) × state.ρu
end

struct GeostrophicForcing{FT} <: Source
  f_coriolis::FT
  u_geostrophic::FT
  v_geostrophic::FT
end
function atmos_source!(s::GeostrophicForcing, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  u = state.ρu / state.ρ
  u_geo = SVector(s.u_geostrophic, s.v_geostrophic, 0)
  source.ρu -= state.ρ * s.f_coriolis * (u - u_geo)
end

"""
  RayleighSponge{FT} <: Source
Rayleigh Damping (Linear Relaxation) for top wall momentum components
Assumes laterally periodic boundary conditions for LES flows. Momentum components
are relaxed to reference values (zero velocities) at the top boundary.
"""
struct RayleighSponge{FT} <: Source
  "Domain maximum height [m]"
  zmax::FT
  "Vertical extent at with sponge starts [m]"
  zsponge::FT
  "Sponge Strength 0 ⩽ c_sponge ⩽ 1"
  c_sponge::FT
  "Ref vel"
  u_ref::SVector{3,FT}
end
function atmos_source!(s::RayleighSponge, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  FT = eltype(state)
  z = aux.orientation.Φ / grav
  coeff = FT(0)
  if z >= s.zsponge
    coeff_top = s.c_sponge * (sinpi(FT(1/2)*(z - s.zsponge)/(s.zmax-s.zsponge)))^FT(4)
    coeff = min(coeff_top, 1.0)
  end
  source.ρu -= (state.ρu - state.ρ * s.u_ref) .* coeff 
end

"""
  Bulk Drag
Testing Parameterisation for DYCOMS with additional momentum reduction term via bulk drag coefficient
"""
struct BulkDrag{FT} <: Source
  "Bulk aerodynamic drag coefficient"
  C_drag::FT
end
function atmos_source!(s::BulkDrag, m::AtmosModel, source::Vars, state::Vars, aux::Vars, t::Real)
  FT = eltype(state)
  u = state.ρu / state.ρ
  z = aux.coord[3]
  if z <= FT(25)
    source.ρu -= state.ρ * s.C_drag * sqrt(sum(abs2.(u))) .* u * (exp(-z/FT(5)))
    source.ρe -= u' * (state.ρ * s.C_drag * sqrt(sum(abs2.(u))) .* u)
  end
end


