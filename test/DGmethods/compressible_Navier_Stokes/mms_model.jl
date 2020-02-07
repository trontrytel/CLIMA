using Unitful; using CLIMA.UnitAnnotations #FIXME
using CLIMA.VariableTemplates

import CLIMA.DGmethods: BalanceLaw, vars_aux, vars_state, vars_gradient,
                        space_unit, time_unit,
                        vars_diffusive, flux_nondiffusive!, flux_diffusive!,
                        source!, wavespeed, boundary_state!,
                        gradvariables!,
                        diffusive!, init_aux!, init_state!,
                        init_ode_state, LocalGeometry


struct MMSModel{dim} <: BalanceLaw
end

sp(T) = U(T, u"m")
ve(T) = U(T, u"m/s")
mf(T) = U(T, u"kg/m^2/s")
ed(T) = U(T, u"J/m^3")

space_unit(m::MMSModel) = u"m"
time_unit(m::MMSModel) = u"s"

vars_aux(::MMSModel,T) = @vars(x1::sp(T), x2::sp(T), x3::sp(T))
vars_state(::MMSModel, T) = @vars(ρ::U(T, u"kg/m^3"), ρu::mf(T), ρv::mf(T), ρw::mf(T), ρe::ed(T))
vars_gradient(::MMSModel, T) = @vars(u::ve(T), v::ve(T), w::ve(T))
vars_diffusive(::MMSModel, T) = @vars(τ11::ed(T), τ22::ed(T), τ33::ed(T), τ12::ed(T), τ13::ed(T), τ23::ed(T))

function flux_nondiffusive!(::MMSModel, flux::Grad, state::Vars,
                            auxstate::Vars, t::Real)
  # preflux
  T = eltype(flux)
  γ = T(γ_exact)
  ρinv = 1 / state.ρ
  u, v, w = ρinv * state.ρu, ρinv * state.ρv, ρinv * state.ρw
  P = (γ-1)*(state.ρe - ρinv * (state.ρu^2 + state.ρv^2 + state.ρw^2) / 2)

  # invisc terms
  flux.ρ  = SVector(state.ρu          , state.ρv          , state.ρw)
  flux.ρu = SVector(u * state.ρu  + P , v * state.ρu      , w * state.ρu)
  flux.ρv = SVector(u * state.ρv      , v * state.ρv + P  , w * state.ρv)
  flux.ρw = SVector(u * state.ρw      , v * state.ρw      , w * state.ρw + P)
  flux.ρe = SVector(u * (state.ρe + P), v * (state.ρe + P), w * (state.ρe + P))
end

function flux_diffusive!(::MMSModel, flux::Grad, state::Vars,
                         diffusive::Vars, auxstate::Vars, t::Real)
  ρinv = 1 / state.ρ
  u, v, w = ρinv * state.ρu, ρinv * state.ρv, ρinv * state.ρw

  # viscous terms
  flux.ρu -= SVector(diffusive.τ11, diffusive.τ12, diffusive.τ13)
  flux.ρv -= SVector(diffusive.τ12, diffusive.τ22, diffusive.τ23)
  flux.ρw -= SVector(diffusive.τ13, diffusive.τ23, diffusive.τ33)

  flux.ρe -= SVector(u * diffusive.τ11 + v * diffusive.τ12 + w * diffusive.τ13,
                     u * diffusive.τ12 + v * diffusive.τ22 + w * diffusive.τ23,
                     u * diffusive.τ13 + v * diffusive.τ23 + w * diffusive.τ33)
end

function gradvariables!(::MMSModel, transformstate::Vars, state::Vars, auxstate::Vars, t::Real)
  ρinv = 1 / state.ρ
  transformstate.u = ρinv * state.ρu
  transformstate.v = ρinv * state.ρv
  transformstate.w = ρinv * state.ρw
end

function diffusive!(m::MMSModel, diffusive::Vars, ∇transform::Grad, state::Vars, auxstate::Vars, t::Real)
  T = eltype(diffusive)
  μ = U(T, u"kg/m/s")(μ_exact)

  dudx, dudy, dudz = ∇transform.u
  dvdx, dvdy, dvdz = ∇transform.v
  dwdx, dwdy, dwdz = ∇transform.w

  # strains
  ϵ11 = dudx
  ϵ22 = dvdy
  ϵ33 = dwdz
  ϵ12 = (dudy + dvdx) / 2
  ϵ13 = (dudz + dwdx) / 2
  ϵ23 = (dvdz + dwdy) / 2

  # deviatoric stresses
  diffusive.τ11 = 2μ * (ϵ11 - (ϵ11 + ϵ22 + ϵ33) / 3)
  diffusive.τ22 = 2μ * (ϵ22 - (ϵ11 + ϵ22 + ϵ33) / 3)
  diffusive.τ33 = 2μ * (ϵ33 - (ϵ11 + ϵ22 + ϵ33) / 3)
  diffusive.τ12 = 2μ * ϵ12
  diffusive.τ13 = 2μ * ϵ13
  diffusive.τ23 = 2μ * ϵ23
end

function source!(m::MMSModel{dim}, source::Vars, state::Vars, aux::Vars, t::Real) where {dim}
  x1,x2,x3 = (aux.x1, aux.x2, aux.x3) ./ space_unit(m)
  source.ρ  = Sρ_g(t, x1, x2, x3, Val(dim))
  source.ρu = SU_g(t, x1, x2, x3, Val(dim))
  source.ρv = SV_g(t, x1, x2, x3, Val(dim))
  source.ρw = SW_g(t, x1, x2, x3, Val(dim))
  source.ρe = SE_g(t, x1, x2, x3, Val(dim))
end

function wavespeed(::MMSModel, nM, state::Vars, aux::Vars, t::Real)
  T = eltype(state)
  γ = T(γ_exact)
  ρinv = 1 / state.ρ
  u, v, w = ρinv * state.ρu, ρinv * state.ρv, ρinv * state.ρw
  P = (γ-1)*(state.ρe - ρinv * (state.ρu^2 + state.ρv^2 + state.ρw^2) / 2)
  return abs(nM[1] * u + nM[2] * v + nM[3] * w) + sqrt(ρinv * γ * P)
end

function boundary_state!(::Rusanov, bl::MMSModel, stateP::Vars, auxP::Vars, nM,
                         stateM::Vars, auxM::Vars, bctype, t, _...)
  init_state!(bl, stateP, auxP, (auxM.x1, auxM.x2, auxM.x3), t)
end

# FIXME: This is probably not right....
boundary_state!(::CentralNumericalFluxGradient, bl::MMSModel, _...) = nothing

function boundary_state!(::CentralNumericalFluxDiffusive, bl::MMSModel,
                         stateP::Vars, diffP::Vars, auxP::Vars, nM,
                         stateM::Vars, diffM::Vars, auxM::Vars, bctype, t, _...)
  init_state!(bl, stateP, auxP, (auxM.x1, auxM.x2, auxM.x3), t)
end

function init_aux!(::MMSModel, aux::Vars, g::LocalGeometry)
  x1,x2,x3 = g.coord
  aux.x1 = x1
  aux.x2 = x2
  aux.x3 = x3
end

function init_state!(bl::MMSModel{dim}, state::Vars, aux::Vars, coord, t) where {dim}
  x1,x2,x3 = coord ./ space_unit(bl)
  state.ρ = ρ_g(t, x1, x2, x3, Val(dim))
  state.ρu = U_g(t, x1, x2, x3, Val(dim))
  state.ρv = V_g(t, x1, x2, x3, Val(dim))
  state.ρw = W_g(t, x1, x2, x3, Val(dim))
  state.ρe = E_g(t, x1, x2, x3, Val(dim))
end
