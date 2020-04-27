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

# rain size distribution parameter
export lambda

# snow intercept parameter
export MP_n_0_sno

# rain fall speed
export c_velo_rai
export terminal_velocity

# rates of conversion between microphysics categories
export conv_q_vap_to_q_liq
export conv_q_vap_to_q_ice
export conv_q_liq_to_q_rai_acnv
export conv_q_liq_to_q_rai_accr
export conv_q_rai_to_q_vap

"""
    c_velo_rai(ρ)

where:
  - `ρ` is the density of air

Returns the proportionality coefficient between terminal velocity of an
individual water drop and the square root of its radius.
"""
function c_velo_rai(ρ::FT) where {FT<:Real}

    # terminal_vel_of_individual_drop = c_velo_rai * drop_radius^(p_velo_rai)
    return sqrt(grav * FT(8/3) / C_drag * (ρ_cloud_liq / ρ - FT(1)))
end

"""
    MP_n_0_sno(q_sno, ρ)

where:
  - `q_sno` is the snow specific humidity
  - `ρ` is the air density

Returns the intercept parameter of the assumed Marshal Palmer distribution of
snow particles.
"""
function MP_n_0_sno(q_sno::FT, ρ::FT) where {FT<:Real}

    return MP_sno_c * (ρ * q_sno)^MP_sno_p
end

"""
    lambda(q_, ρ, c_mass, p_mass)

where:
  - `ρ` is the density of air
  - `q_` is specific humidity of rain, ice or snow
  - `MP_n_0` is the intercept parameter of the assumed size distribution
  - `c_mass` is the coefficient in the mass(radius) relationship
  - `p_mass` is the exponent in the mass(radius) relationship

Returns the rate parameter of the assumed size distribution of
particles (rain drops, ice crystals, snow crystals).
"""
function lambda(q_::FT, ρ::FT, MP_n_0::FT, c_mass::FT, p_mass::FT) where {FT<:Real}

    return (c_mass * MP_n_0 * gamma(p_mass + FT(1)) / ρ / q_rai)^FT(1/(p_mass + 1))
end

"""
    terminal_velocity(q_rai, ρ, MP_n_0, c_velo, p_velo, c_mass, p_mass)

where:
  - `q_` is the specific humidity of rain or snow
  - `ρ`  is the density of air
  - `MP_n_0` is the intercept parameter of the assumed size distribution for rain and snow
  - `c_velo, p_velo` is the velocity(radius) coefficients
  - `c_mass, p_mass` is the mass(radius) coefficients

Returns the mass weighted average terminal velocity assuming
Marshall Palmer 1948 distribution of rain drops and snow crystals.
"""
function terminal_velocity(q_::FT, ρ::FT, MP_n_0::FT, c_velo::FT, p_velo::FT,
                           c_mass::FT, p_mass::FT) where {FT <: Real}

    fall_w = FT(0)

    if q_ > FT(0)

        λ_ = lambda(q_, ρ, MP_n_0, c_mass, p_mass)

        fall_w = c_velo * λ_^(-p_velo) *
                 gamma(p_velo + p_mass + FT(1)) / gamma(p_mass + FT(1))
    end
    return fall_w
end

"""
    conv_q_vap_to_q_liq(q_sat, q)

where:
- `q_sat` - PhasePartition at equilibrium
- `q`     - current PhasePartition

Returns the q_liq tendency due to condensation/evaporation.
The tendency is obtained assuming a relaxation to equilibrium with
constant timescale.
"""
function conv_q_vap_to_q_liq(q_sat::PhasePartition{FT},
                             q::PhasePartition{FT}) where {FT<:Real}
    return (q_sat.liq - q.liq) / τ_cond_evap_liq
end

"""
    conv_q_vap_to_q_ice(q_sat, q)

where:
- `q_sat` - PhasePartition at equilibrium
- `q`     - current PhasePartition

Returns the q_ice tendency due to sublimation/resublimation.
The tendency is obtained assuming a relaxation to equilibrium with
constant timescale.
"""
function conv_q_vap_to_q_ice(q_sat::PhasePartition{FT},
                             q::PhasePartition{FT}) where {FT<:Real}
    return (q_sat.ice - q.ice) / τ_cond_evap_ice
end

"""
    conv_q_liq_to_q_rai_acnv(q_liq)

where:
- `q_liq` - is the liquid water specific humidity

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion) parametrized following Kessler 1995.
"""
function conv_q_liq_to_q_rai_acnv(q_liq::FT) where {FT <: Real}

    return max(FT(0), q_liq - q_liq_threshold) / τ_acnv
end


"""
    accretion(q_susp, q_fsll, ρ)

where:
- `q_susp` is the suspended water specific humidity (cloud or ice)
- `q_fall` is the precipitating water specific humidity (rain or snow)
- `ρ` is the density of air
- `MP_n_0`
- α_geom
- c_area, p_area
- c_mass, p_mass
- c_velo, p_velo

Returns the q_rai tendency due to collisions between cloud droplets
and rain drops (accretion) parametrized following Kessler 1995.
"""
function accretion(q_susp::FT, q_fall::FT, ρ::FT,
                   E_col::FT, MP_n_0::FT, α_geom::FT,
                   c_area::FT, p_area::FT, c_velo::FT, p_velo::FT,
                   c_mass::FT, p_mass::FT) where {FT<:Real}

    accr_rate = FT(0)

    if (q_rai > FT(0) && q_liq > FT(0))

        λ_ =  lambda(q_fall, ρ, MP_n_0, c_mass, p_mass)

        accr_rate = MP_n_0 * c_area * c_velo(ρ) * q_susp * E_col * α_geom *
                    gamma(p_area + p_velo + FT(1)) *
                    λ_^(-p_area - p_velo -FT(1))
    end
    return accr_rate
end

"""
    conv_q_rai_to_q_vap(q_rai, q, T, p, ρ)

where:
 - q_rai - rain water specific humidity
 - q - current PhasePartition
 - T - temperature
 - p - pressure
 - ρ - air density

Returns the q_rai tendency due to rain evaporation. Parameterized following
Smolarkiewicz and Grabowski 1996.
"""
function conv_q_rai_to_q_vap(q_rai::FT, q::PhasePartition{FT},
                             T::FT, p::FT, ρ::FT) where {FT<:Real}
    evap_rate = FT(0)

    if q_rai > FT(0)

        qv_sat = q_vap_saturation(T, ρ, Liquid())
        q_v = q.tot - q.liq - q.ice
        S = q_v/qv_sat - FT(1)

        L = latent_heat_vapor(T)
        p_vs = saturation_vapor_pressure(T, Liquid())
        G::FT = FT(1) / (
                  L / K_therm / T * (L / R_v / T - FT(1))
                  + R_v * T / D_vapor / p_vs
                )

        λ_rai =  lambda(q_rai, ρ, MP_n_0_rai, c_mass_rai, p_mass_rai)

        evap_rate = c_mass_rai * p_mass_rai * S * G * MP_n_0_rai / ρ / ρ_cloud_liq *
                    (a_vent * λ_rai^(-p_mass_rai + FT(1)) * gamma(p_mass_rai - FT(1))
                     +
                     b_vent * N_Sc^FT(1/3) * (FT(2) * c_velo_rai(ρ) / ν_air *
                     λ_rai^FT(-2*p_mass_rai - p_velo_rai + 1))^FT(1/2) *
                     gamma(FT(2*p_mass_rai + p_velo_rai - 1)/FT(2))
                    )
    end
    return evap_rate
end

"""
    conv_q_ice_to_q_sno(q, T, p, ρ)

where:
 - q - current PhasePartition
 - T - temperature
 - p - pressure
 - ρ - air density

Returns the q_sno tendency due to autoconversion from ice.
Parameterized following Harrington et al 1996 and Kaul et al 2015
"""
function conv_q_ice_to_q_sno(q::PhasePartition{FT},
                             T::FT, p::FT, ρ::FT) where {FT<:Real}
    acnv_rate = FT(0)

    if q.ice > FT(0)

        qv_sat = q_vap_saturation(T, ρ, Ice())
        q_v = q.tot - q.liq - q.ice
        S = q_v/qv_sat - FT(1)

        L = latent_heat_sublim(T)
        p_vs = saturation_vapor_pressure(T, Ice())
        G::FT = FT(1) / (
                  L / K_therm / T * (L / R_v / T - FT(1))
                  + R_v * T / D_vapor / p_vs
                )

        λ_ice =  lambda(q.ice, ρ, MP_n_0_ice, c_mass_ice, p_mass_ice)

        acnv_rate = FT(4) * π * S * G * MP_n_0_ice / ρ *
                    exp(-λ_ice * r_ice_snow) *
                    (
                      r_ice_snow^FT(2) / p_mass_ice
                      +
                      (r_ice_snow * λ_ice + FT(1)) / λ_ice^FT(2)
                    )

    end
    return acnv_rate
end





end #module Microphysics.jl
