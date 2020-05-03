"""
    one-moment bulk Microphysics scheme

Microphysics parameterization based on the ideas of Kessler_1995:
  - condensation/evaporation as relaxation to equilibrium
  - autoconversion
  - accretion
  - rain evaporation
  - rain terminal velocity
"""
module Microphysics

using SpecialFunctions

using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters: ρ_cloud_liq, R_v, grav
using ..MicrophysicsParameters
using ..Parameters

export ζ_rai
export n0_sno
export lambda

export supersaturation
export G_func

export terminal_velocity

export conv_q_vap_to_q_liq
export conv_q_vap_to_q_ice

export conv_q_liq_to_q_rai_acnv
export conv_q_ice_to_q_sno_acnv

export accretion
export accretion_rain_sink
export accretion_snow_rain

export evaporation_sublimation
export snow_melt

"""
    ζ_rai(ρ)

 - `ρ` air density

Returns the proportionality coefficient between terminal velocity of an
individual water drop and the square root of its radius.
"""
function ζ_rai(ρ::FT) where {FT<:Real}

    return sqrt(grav * FT(8/3) / C_drag * (ρ_cloud_liq / ρ - FT(1)))
end

"""
    n0_sno(q_sno, ρ)

  - `q_sno` -  snow specific humidity
  - `ρ` - air density

Returns the intercept parameter of the assumed Marshal Palmer distribution of
snow particles.
"""
function n0_sno(q_sno::FT, ρ::FT) where {FT<:Real}

    return μ_sno * (ρ * q_sno)^ν_sno
end

"""
    lambda(q, ρ, n0, α, β)

  - `q` - specific humidity of rain, ice or snow
  - `ρ` - air density
  - `n0` - intercept parameter of the assumed size distribution
  - `α`, `β` - mass(radius) parameters

Returns the rate parameter of the assumed size distribution of
particles (rain drops, ice crystals, snow crystals).
"""
function lambda(q::FT, ρ::FT, n0::FT, α::FT, β::FT) where {FT<:Real}

    return (α * n0 * gamma(β + FT(1)) / ρ / q)^FT(1/(β + 1))
end

"""
    supersaturation(q, ρ, Liquid())
    supersaturation(q, ρ, Ice())

 - `q` - phase partition
 - `ρ` - air density,
 - `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.

Returns supersaturation (qv/qv_sat -1) over water or ice.
"""
function supersaturation(q::PhasePartition{FT}, ρ::FT,
                         ::Liquid) where {FT<:Real}

    q_sat = q_vap_saturation(T, ρ, Liquid())
    q_vap = q.tot - q.liq - q.ice

    return q_vap/q_sat - FT(1)
end
function supersaturation(q::PhasePartition{FT}, ρ::FT, ::Ice) where {FT<:Real}

    q_sat = q_vap_saturation(T, ρ, Ice())
    q_vap = q.tot - q.liq - q.ice

    return q_vap/q_sat - FT(1)
end

"""
    G_func(T, K_therm, R_v, D_vapor, Liquid())
    G_func(T, K_therm, R_v, D_vapor, Ice())

 - `T` - air temperature
 - `K_therm` - thermal conductivity of air
 - `R_v` - water vapor gas constant
 - `D_vapor` - water vapor diffusivity
 - `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.

Utility function combining thermal conductivity and vapor diffusivity effects.
"""
function G_func(T::FT, K_therm::FT, R_v::FT, D_vapor::FT,
                ::Liquid) where {FT<:Real}

    L = latent_heat_vapor(T)
    p_vs = saturation_vapor_pressure(T, Liquid())

    return FT(1) / (
              L / K_therm / T * (L / R_v / T - FT(1))
              + R_v * T / D_vapor / p_vs
           )
end
function G_func(T::FT, K_therm::FT, R_v::FT, D_vapor::FT,
                ::Ice) where {FT<:Real}

    L = latent_heat_sublim(T)
    p_vs = saturation_vapor_pressure(T, Ice())

    return FT(1) / (
              L / K_therm / T * (L / R_v / T - FT(1))
              + R_v * T / D_vapor / p_vs
           )
end

"""
    terminal_velocity(q, λ, β, ζ, η)

  - `q` - rain or snow specific humidity
  - `λ` - size distribution parameter
  - `β` - mass(radius) coefficient
  - `ζ`, `η` - velocity(radius) coefficients

Returns the mass weighted average terminal velocity assuming
Marshall Palmer 1948 distribution of rain drops and snow crystals.
"""
function terminal_velocity(q::FT, λ::FT, β::FT, ζ::FT,
                           η::FT) where {FT <: Real}
    fall_w = FT(0)

    if q > FT(0)
        fall_w = ζ * λ^(-η) * gamma(η + β + FT(1)) / gamma(β + FT(1))
    end

    return fall_w
end

"""
    conv_q_vap_to_q_liq_ice(q_sat, q, Liquid())
    conv_q_vap_to_q_liq_ice(q_sat, q, Ice())

- `q_sat` - PhasePartition at equilibrium
- `q`     - current PhasePartition
- `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.

Returns the q_liq tendency due to condensation/evaporation
or q_ice tendency due to sublimation/resublimation.

The tendency is obtained assuming a relaxation to equilibrium with
constant timescale.
"""
function conv_q_vap_to_q_liq_ice(q_sat::PhasePartition{FT},
                             q::PhasePartition{FT}, ::Liquid) where {FT<:Real}
    return (q_sat.liq - q.liq) / τ_cond_evap
end
function conv_q_vap_to_q_liq_ice(q_sat::PhasePartition{FT},
                             q::PhasePartition{FT}, ::Ice) where {FT<:Real}
    return (q_sat.ice - q.ice) / τ_sub_resub
end

"""
    conv_q_liq_to_q_rai(q_liq)

- `q_liq` - liquid water specific humidity

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion) parametrized following Kessler 1995.
"""
function conv_q_liq_to_q_rai(q_liq::FT) where {FT <: Real}

    return max(FT(0), q_liq - q_liq_threshold) / τ_acnv
end

"""
    conv_q_ice_to_q_sno(q_ice, ρ, S, G, r_ice_sno, λ_ice, n0_ice, β_ice)

 - `q_ice` - snow specific humidity
 - `ρ` - air density
 - `S` - supersaturation
 - `G` - utility function `G_func`
 - `r_ice_sno` - threshold radius between ice and snow
 - `λ_ice`m `n0_ice` - ice size distribution parameters
 - `β_ice` - ice mass(radius) parameter

Returns the q_sno tendency due to autoconversion from ice.
Parameterized following Harrington et al 1996 and Kaul et al 2015
"""
function conv_q_ice_to_q_sno_acnv(q_ice::FT, ρ::FT, S::FT, G::FT,
                                  r_ice_sno::FT, λ_ice::FT, n0_ice::FT,
                                  β_ice::FT) where {FT<:Real}
    acnv_rate = FT(0)

    if q_ice > FT(0)

        acnv_rate = FT(4) * π * S * G * n0_ice / ρ *
                    exp(-λ_ice * r_ice_snow) *
                    (
                      r_ice_snow^FT(2) / β_ice
                      +
                      (r_ice_snow * λ_ice + FT(1)) / λ_ice^FT(2)
                    )

    end
    return acnv_rate
end

"""
    accretion(q_susp, q_fall, ρ, E_sf, n0_fall,
              α_fall, β_fall, γ_fall, δ_fall, ζ_fall, η_fall)

 - `q_susp` is the suspended water specific humidity (cloud or ice)
 - `q_fall` is the precipitating water specific humidity (rain or snow)
 - `E_sf` - is the collision efficiency
 - `n0_fall`, `λ_fall` - is the size distribution coefficient
 - `γ_fall`, `δ_fall`, `ζ_fall`, `η_fall` are the size relationships coefficients

Returns the sink to suspended water (cloud water or cloud ice) due to collisions
with falling water (rain or snow).
"""
function accretion(q_susp::FT, q_fall::FT,
                   E_sf::FT, n0_fall::FT, λ_fall::FT,
                   γ_fall::FT, δ_fall::FT,
                   ζ_fall::FT, η_fall::FT) where {FT<:Real}

    accr_rate = FT(0)

    if (q_susp > FT(0) && q_fall > FT(0))

        accr_rate = q_susp * E_sf * n0_fall * γ_fall * ζ_fall *
                    gamma(δ_fall + η_fall + FT(1)) /
                    λ_fall^(δ_fall + η_fall + FT(1))
    end
    return accr_rate
end

"""
    accretion_rain_sink(q_ice, q_rai, ρ, E_ir, n0_ice, n0_rai, α_ice, β_ice,
                        α_rai, β_rai, γ_rai, δ_rai, ζ_rai, η_rai)

 - `q_ice` is the suspended water specific humidity (cloud or ice)
 - `q_rai` is the precipitating water specific humidity (rain or snow)
 - `ρ` is the density of air
 - `E_ir` - is the collision efficiency for rain - ice collisions
 - `n0_ice` and `n0_rai` - are the size distribution coefficients for cloud ice and rain
 - `α_rai`, `β_rai`, `γ_rai`, `δ_rai`, `ζ_rai`, `η_rai` are the size relationships coefficients for rain
 - `α_ice`, `β_ice` are the size relationships coefficients for ice

Returns the sink of rain water (partial source of snow) due to collisions
with cloud ice.
"""
function accretion_rain_sink(q_ice::FT, q_rai::FT, E_ir::FT,
                             n0_ice::FT, n0_rai::FT, λ_ice::FT, λ_rai::FT,
                             α_rai::FT, β_rai::FT
                             γ_rai::FT, δ_rai::FT,
                             ζ_rai::FT, η_rai::FT) where {FT<:Real}
tutaj
    accr_rate = FT(0)
    if (q_ice > FT(0) && q_rai > FT(0))

        accr_rate = E_ir * n0_rai * n0_ice * α_rai * γ_rai * ζ_rai / λ_ice
                    gamma(β_rai + δ_rai + η_rai + FT(1)) *
                    λ_^(-β_rai - δ_rai - η_rai - FT(1))
    end
    return accr_rate
end

"""
    accretion_snow_rain(q_i, q_j, ρ, E_ij, n0_i, n0_j,
                        α_i, β_i, ζ_i, η_i,
                        α_j, β_j, ζ_j, η_j)

where:
- `q_i` - q_snow if T > T_freezing or q_rain if T < T_freezing
- `q_j` - q_rain if T > T_freezing or q_snow if T < T_freezing
- `ρ` is the density of air
- `E_ij` - is the collision efficiency for rain - snow collisions
- `n0_i` and `n0_j` - are the size distribution coefficients for rain and snow
- `α_i`, `β_i`, `ζ_i`, `η_i` are the size relationships coefficients for snow if T > T_freezing or rain if T < T_freezing
- `α_j`, `β_j`, `ζ_j`, `η_j` are the size relationships coefficients for rain if T > T_freezing or snow if T < T_freezing

Returns the accretion rate between rain and snow.
"""
function accretion_snow_rain(q_i::FT, q_j::FT, ρ::FT,
                             E_ij::FT, n0_i::FT, n0_j::FT,
                             α_i::FT, β_i::FT, α_j::FT, β_j::FT,
                             ζ_i::FT, η_i::FT, ζ_j::FT, η_j::FT
                            ) where {FT<:Real}

    accr_rate = FT(0)
    if (q_i > FT(0) && q_j > FT(0))

        λ_i =  lambda(q_i, ρ, n0_i, α_i, β_i)
        λ_j =  lambda(q_j, ρ, n0_j, α_j, β_j)

        v_ti = terminal_velocity(q_i, ρ, n0_i, α_i, β_i, ζ_i, η_i)
        v_tj = terminal_velocity(q_j, ρ, n0_j, α_j, β_j, ζ_j, η_j)

        accr_rate = π * n0_i * n0_j * α_j * E_ij * abs(v_ti - v_tj) *
                    (
                      FT(2) * gamma(β_j + FT(1)) / λ_i^FT(3) / λ_j^(β_j + FT(1)) +
                      FT(2) * gamma(β_j + FT(2)) / λ_i^FT(2) / λ_j^(β_j + FT(2)) +
                      gamma(β_j + FT(3)) / λ_i /λ_j^(β_j + FT(3))
                    )
    end
    return accr_rate
end

"""
    evaporation_sublimation(q, ρ, n0, S, G, λ, a_vent, b_vent, ν_air, D_vapor,
                            λ, ζ, η)
where:
 - q - rain or snow water specific humidity
 - ρ - air density
 - n0, λ - size distribution parameters
 - S - supersaturation
 - G - function combining the effects of thermal conductivity and vapor diffusivity
       (see snow autoconversion docs for details)
 - a_vent, b_vent - ventilation factor coefficients
 - ν_air - kinematic viscosity of air
 - D_vapor - diffusity of water vapor
 - ζ, η - terminal velocity parameters

Returns the tendency due to rain evaporation or snow sublimation.
"""
function evaporation_sublimation(q::FT, ρ::FT, n0::FT, S::FT, G::FT, λ::FT,
                                 a_vent::FT, b_vent::FT, ν_air::FT, D_vapor::FT,
                                 λ::FT, ζ::FT, η::FT) where {FT<:Real}
    evap_subl_rate = FT(0)

    if q  > FT(0)

        evap_subl_rate = 4 * π * n0 / ρ * S * G / λ^FT(2) * (
          a_vent +
          b_vent * (ν_air / D_vapor)^FT(1/3) * (FT(2) * ζ / ν_air)^FT(1/2) *
            gamma((η + FT(5)) / FT(2)) + λ^(-(λ + FT(1))/FT(2))
        )


    end
    return evap_subl_rate
end

"""
    snow_melt(q_sno, ρ, n0_sno, λ_sno, a_vent_sno, b_vent_sno,
              ν_air, D_vapor, K_therm, L, T_freeze, λ_sno, ζ_sno, η_sno)
where:
 - q_sno - snow water specific humidity
 - ρ - air density
 - n0_sno, λ_sno - size distribution parameters
 - a_vent_sno, b_vent_sno - ventilation factor coefficients
 - ν_air - kinematic viscosity of air
 - D_vapor - diffusity of water vapor
 - K_them - thermal conductivity of air
 - L - latent heat of freezing
 - T_freeze - freezing temperature
 - ζ_sno, η_sno - terminal velocity parameters

Returns the tendency due to rain evaporation or snow sublimation.
"""
function snow_melt(q_sno::FT, ρ::FT, n0_sno::FT, λ_sno::FT,
                   a_vent_sno::FT, b_vent_sno::FT, ν_air::FT, D_vapor::FT,
                   K_therm::FT, L::FT, T_freeze::FT,
                   λ_sno::FT, ζ_sno::FT, η_sno::FT) where {FT<:Real}

    snow_melt_rate = FT(0)

    if q_sno  > FT(0)

        snow_melt_rate = 4 * π * n0_sno / ρ * K_therm / L * (T - T_freeze) / λ^FT(2) * (
          a_vent +
          b_vent * (ν_air / D_vapor)^FT(1/3) * (FT(2) * ζ / ν_air)^FT(1/2) *
            gamma((η + FT(5)) / FT(2)) + λ^(-(λ + FT(1))/FT(2))
        )


    end
    return snow_melt_rate
end

end #module Microphysics.jl
