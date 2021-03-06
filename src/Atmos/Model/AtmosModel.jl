module Atmos

export AtmosModel,
    AtmosAcousticLinearModel, AtmosAcousticGravityLinearModel, RemainderModel

using LinearAlgebra, StaticArrays
using ..ConfigTypes
using ..VariableTemplates
using ..Parameters
using ..MoistThermodynamics
using ..PlanetParameters
import ..MoistThermodynamics: internal_energy
using ..SubgridScaleParameters
using ..MPIStateArrays: MPIStateArray
using ..Mesh.Grids: VerticalDirection, HorizontalDirection, min_node_distance

import CLIMA.DGmethods:
    BalanceLaw,
    vars_aux,
    vars_state,
    vars_gradient,
    vars_gradient_laplacian,
    vars_diffusive,
    vars_hyperdiffusive,
    flux_nondiffusive!,
    flux_diffusive!,
    source!,
    wavespeed,
    boundary_state!,
    gradvariables!,
    diffusive!,
    hyperdiffusive!,
    init_aux!,
    init_state!,
    update_aux!,
    LocalGeometry,
    lengthscale,
    resolutionmetric,
    DGModel,
    nodal_update_aux!,
    num_state,
    num_integrals,
    vars_integrals,
    vars_reverse_integrals,
    indefinite_stack_integral!,
    reverse_indefinite_stack_integral!,
    integral_load_aux!,
    integral_set_aux!,
    reverse_integral_load_aux!,
    reverse_integral_set_aux!
import ..DGmethods.NumericalFluxes:
    boundary_state!,
    boundary_flux_diffusive!,
    normal_boundary_flux_diffusive!,
    NumericalFluxNonDiffusive,
    NumericalFluxGradient,
    NumericalFluxDiffusive,
    CentralHyperDiffusiveFlux,
    CentralDivPenalty

import ..Courant: advective_courant, nondiffusive_courant, diffusive_courant

"""
    AtmosModel <: BalanceLaw

A `BalanceLaw` for atmosphere modeling.

# Usage

    AtmosModel(orientation, ref_state, turbulence, hyperdiffusion, moisture, radiation, source,
               boundarycondition, init_state)

"""
struct AtmosModel{FT, PS, O, RS, T, HD, M, P, R, S, BC, IS, DC} <: BalanceLaw
    param_set::PS
    orientation::O
    ref_state::RS
    turbulence::T
    hyperdiffusion::HD
    moisture::M
    precipitation::P
    radiation::R
    source::S
    # TODO: Probably want to have different bc for state and diffusion...
    boundarycondition::BC
    init_state::IS
    data_config::DC
end

function AtmosModel{FT}(
    ::Type{AtmosLESConfigType};
    orientation::O = FlatOrientation(),
    ref_state::RS = HydrostaticState(
        LinearTemperatureProfile(FT(200), FT(280), FT(grav) / FT(cp_d)),
        FT(0),
    ),
    turbulence::T = SmagorinskyLilly{FT}(0.21),
    hyperdiffusion::HD = NoHyperDiffusion(),
    moisture::M = EquilMoist{FT}(),
    precipitation::P = NoPrecipitation(),
    radiation::R = NoRadiation(),
    source::S = (Gravity(), Coriolis(), GeostrophicForcing{FT}(7.62e-5, 0, 0)),
    # TODO: Probably want to have different bc for state and diffusion...
    boundarycondition::BC = AtmosBC(),
    init_state::IS = nothing,
    data_config::DC = nothing,
    param_set::PS = nothing,
) where {FT <: AbstractFloat, O, RS, T, HD, M, P, R, S, BC, IS, DC, PS}
    @assert param_set ≠ nothing
    @assert init_state ≠ nothing

    atmos = (
        param_set,
        orientation,
        ref_state,
        turbulence,
        hyperdiffusion,
        moisture,
        precipitation,
        radiation,
        source,
        boundarycondition,
        init_state,
        data_config,
    )

    return AtmosModel{FT, typeof.(atmos)...}(atmos...)
end
function AtmosModel{FT}(
    ::Type{AtmosGCMConfigType};
    orientation::O = SphericalOrientation(),
    ref_state::RS = HydrostaticState(
        LinearTemperatureProfile(FT(200), FT(280), FT(grav) / FT(cp_d)),
        FT(0),
    ),
    turbulence::T = SmagorinskyLilly{FT}(0.21),
    hyperdiffusion::HD = NoHyperDiffusion(),
    moisture::M = EquilMoist{FT}(),
    precipitation::P = NoPrecipitation(),
    radiation::R = NoRadiation(),
    source::S = (Gravity(), Coriolis()),
    boundarycondition::BC = AtmosBC(),
    init_state::IS = nothing,
    data_config::DC = nothing,
    param_set::PS = nothing,
) where {FT <: AbstractFloat, O, RS, T, HD, M, P, R, S, BC, IS, DC, PS}
    @assert param_set ≠ nothing
    @assert init_state ≠ nothing
    atmos = (
        param_set,
        orientation,
        ref_state,
        turbulence,
        hyperdiffusion,
        moisture,
        precipitation,
        radiation,
        source,
        boundarycondition,
        init_state,
        data_config,
    )

    return AtmosModel{FT, typeof.(atmos)...}(atmos...)
end

function vars_state(m::AtmosModel, FT)
    @vars begin
        ρ::FT
        ρu::SVector{3, FT}
        ρe::FT
        turbulence::vars_state(m.turbulence, FT)
        hyperdiffusion::vars_state(m.hyperdiffusion, FT)
        moisture::vars_state(m.moisture, FT)
        radiation::vars_state(m.radiation, FT)
    end
end
function vars_gradient(m::AtmosModel, FT)
    @vars begin
        u::SVector{3, FT}
        h_tot::FT
        turbulence::vars_gradient(m.turbulence, FT)
        hyperdiffusion::vars_gradient(m.hyperdiffusion, FT)
        moisture::vars_gradient(m.moisture, FT)
    end
end
function vars_diffusive(m::AtmosModel, FT)
    @vars begin
        ∇h_tot::SVector{3, FT}
        turbulence::vars_diffusive(m.turbulence, FT)
        hyperdiffusion::vars_diffusive(m.hyperdiffusion, FT)
        moisture::vars_diffusive(m.moisture, FT)
    end
end
function vars_gradient_laplacian(m::AtmosModel, FT)
    @vars begin
        hyperdiffusion::vars_gradient_laplacian(m.hyperdiffusion, FT)
    end
end
function vars_hyperdiffusive(m::AtmosModel, FT)
    @vars begin
        hyperdiffusion::vars_hyperdiffusive(m.hyperdiffusion, FT)
    end
end
function vars_aux(m::AtmosModel, FT)
    @vars begin
        ∫dz::vars_integrals(m, FT)
        ∫dnz::vars_reverse_integrals(m, FT)
        coord::SVector{3, FT}
        orientation::vars_aux(m.orientation, FT)
        ref_state::vars_aux(m.ref_state, FT)
        turbulence::vars_aux(m.turbulence, FT)
        hyperdiffusion::vars_aux(m.hyperdiffusion, FT)
        moisture::vars_aux(m.moisture, FT)
        radiation::vars_aux(m.radiation, FT)
    end
end
function vars_integrals(m::AtmosModel, FT)
    @vars begin
        radiation::vars_integrals(m.radiation, FT)
    end
end
function vars_reverse_integrals(m::AtmosModel, FT)
    @vars begin
        radiation::vars_reverse_integrals(m.radiation, FT)
    end
end

include("orientation.jl")
include("ref_state.jl")
include("turbulence.jl")
include("hyperdiffusion.jl")
include("moisture.jl")
include("precipitation.jl")
include("radiation.jl")
include("source.jl")
include("boundaryconditions.jl")
include("linear.jl")
include("remainder.jl")
include("courant.jl")

"""
    flux_nondiffusive!(m::AtmosModel, flux::Grad, state::Vars, aux::Vars,
                       t::Real)

Computes flux non-diffusive flux portion of `F` in:

```
∂Y
-- = - ∇ • (F_{adv} + F_{press} + F_{nondiff} + F_{diff}) + S(Y)
∂t
```
Where

 - `F_{adv}`      Advective flux             ; see [`flux_advective!`]@ref()
 - `F_{press}`    Pressure flux              ; see [`flux_pressure!`]@ref()
 - `F_{diff}`     Fluxes that state gradients; see [`flux_diffusive!`]@ref()
"""
@inline function flux_nondiffusive!(
    m::AtmosModel,
    flux::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ρ = state.ρ
    ρinv = 1 / ρ
    ρu = state.ρu
    u = ρinv * ρu

    # advective terms
    flux.ρ = ρ * u
    flux.ρu = ρ * u .* u'
    flux.ρe = u * state.ρe

    # pressure terms
    p = pressure(m, m.moisture, state, aux)
    if m.ref_state isa HydrostaticState
        flux.ρu += (p - aux.ref_state.p) * I
    else
        flux.ρu += p * I
    end
    flux.ρe += u * p
    flux_radiation!(m.radiation, m, flux, state, aux, t)
    flux_moisture!(m.moisture, m, flux, state, aux, t)
end

function gradvariables!(
    atmos::AtmosModel,
    transform::Vars,
    state::Vars,
    aux::Vars,
    t::Real,
)
    ρinv = 1 / state.ρ
    transform.u = ρinv * state.ρu
    transform.h_tot = total_specific_enthalpy(atmos, atmos.moisture, state, aux)

    gradvariables!(atmos.moisture, transform, state, aux, t)
    gradvariables!(atmos.turbulence, transform, state, aux, t)
    gradvariables!(atmos.hyperdiffusion, transform, state, aux, t)
end

function diffusive!(
    atmos::AtmosModel,
    diffusive::Vars,
    ∇transform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    diffusive.∇h_tot = ∇transform.h_tot

    # diffusion terms required for SGS turbulence computations
    diffusive!(
        atmos.turbulence,
        atmos.orientation,
        diffusive,
        ∇transform,
        state,
        aux,
        t,
    )
    # diffusivity of moisture components
    diffusive!(atmos.moisture, diffusive, ∇transform, state, aux, t)
end

function hyperdiffusive!(
    atmos::AtmosModel,
    hyperdiffusive::Vars,
    hypertransform::Grad,
    state::Vars,
    aux::Vars,
    t::Real,
)
    hyperdiffusive!(
        atmos.hyperdiffusion,
        hyperdiffusive,
        hypertransform,
        state,
        aux,
        t,
    )
end

@inline function flux_diffusive!(
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    diffusive::Vars,
    hyperdiffusive::Vars,
    aux::Vars,
    t::Real,
)
    ν, τ = turbulence_tensors(atmos.turbulence, state, diffusive, aux, t)
    D_t = (ν isa Real ? ν : diag(ν)) * inv_Pr_turb
    d_h_tot = -D_t .* diffusive.∇h_tot
    flux_diffusive!(atmos, flux, state, τ, d_h_tot)
    flux_diffusive!(atmos.moisture, flux, state, diffusive, aux, t, D_t)
    flux_diffusive!(
        atmos.hyperdiffusion,
        flux,
        state,
        diffusive,
        hyperdiffusive,
        aux,
        t,
    )
end

#TODO: Consider whether to not pass ρ and ρu (not state), foc BCs reasons
@inline function flux_diffusive!(
    atmos::AtmosModel,
    flux::Grad,
    state::Vars,
    τ,
    d_h_tot,
)
    flux.ρu += τ * state.ρ
    flux.ρe += τ * state.ρu
    flux.ρe += d_h_tot * state.ρ
end

@inline function wavespeed(m::AtmosModel, nM, state::Vars, aux::Vars, t::Real)
    ρinv = 1 / state.ρ
    u = ρinv * state.ρu
    return abs(dot(nM, u)) + soundspeed(m, m.moisture, state, aux)
end


function update_aux!(dg::DGModel, m::AtmosModel, Q::MPIStateArray, t::Real)
    FT = eltype(Q)
    auxstate = dg.auxstate

    if num_integrals(m, FT) > 0
        indefinite_stack_integral!(dg, m, Q, auxstate, t)
        reverse_indefinite_stack_integral!(dg, m, Q, auxstate, t)
    end

    nodal_update_aux!(atmos_nodal_update_aux!, dg, m, Q, t)

    return true
end

function atmos_nodal_update_aux!(m::AtmosModel, state::Vars, aux::Vars, t::Real)
    atmos_nodal_update_aux!(m.moisture, m, state, aux, t)
    atmos_nodal_update_aux!(m.radiation, m, state, aux, t)
    atmos_nodal_update_aux!(m.turbulence, m, state, aux, t)
end

function integral_load_aux!(m::AtmosModel, integ::Vars, state::Vars, aux::Vars)
    integral_load_aux!(m.radiation, integ, state, aux)
end

function integral_set_aux!(m::AtmosModel, aux::Vars, integ::Vars)
    integral_set_aux!(m.radiation, aux, integ)
end

function reverse_integral_load_aux!(
    m::AtmosModel,
    integ::Vars,
    state::Vars,
    aux::Vars,
)
    reverse_integral_load_aux!(m.radiation, integ, state, aux)
end

function reverse_integral_set_aux!(m::AtmosModel, aux::Vars, integ::Vars)
    reverse_integral_set_aux!(m.radiation, aux, integ)
end

# TODO: figure out a nice way to handle this
function init_aux!(m::AtmosModel, aux::Vars, geom::LocalGeometry)
    aux.coord = geom.coord
    atmos_init_aux!(m.orientation, m, aux, geom)
    atmos_init_aux!(m.ref_state, m, aux, geom)
    atmos_init_aux!(m.turbulence, m, aux, geom)
    atmos_init_aux!(m.hyperdiffusion, m, aux, geom)
end

"""
    source!(m::AtmosModel, source::Vars, state::Vars, diffusive::Vars, aux::Vars, t::Real)
Computes flux `S(Y)` in:
```
∂Y
-- = - ∇ • F + S(Y)
∂t
```
"""
function source!(
    m::AtmosModel,
    source::Vars,
    state::Vars,
    diffusive::Vars,
    aux::Vars,
    t::Real,
)
    atmos_source!(m.source, m, source, state, diffusive, aux, t)
end

function init_state!(m::AtmosModel, state::Vars, aux::Vars, coords, t, args...)
    m.init_state(m, state, aux, coords, t, args...)
end
end # module
