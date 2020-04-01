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

using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters: ρ_cloud_liq, R_v, grav
using ..MicrophysicsParameters
using ..Parameters

# rain fall speed
export c_velo
export terminal_velocity

# rain size distribution parameter
export lambda

# rates of conversion between microphysics categories
export conv_q_vap_to_q_liq
export conv_q_vap_to_q_ice
export conv_q_liq_to_q_rai_acnv
export conv_q_liq_to_q_rai_accr
export conv_q_rai_to_q_vap

"""
    c_velo(ρ)

where:
  - `ρ` is the density of air

Returns the proportionality coefficient between terminal velocity of an
individual water drop and the square root of its radius.
"""
function c_velo(ρ::FT) where {FT<:Real}

    # terminal_vel_of_individual_drop = c_velo * drop_radius^(1/2)
    return sqrt(grav * FT(8/3) / C_drag * (ρ_cloud_liq / ρ - FT(1)))
end

"""
    lambda(q_rai, ρ)

where:
  - `ρ`     - the density of air
  - `q_rai` - rain water specific humidity
Returns the coefficient in the exponent of the assumed size distribution of
rain drops.
"""
function lambda(q_rai::FT, ρ::FT) where {FT<:Real}

    return
      (c_mass * MP_n_0 * gamma(p_mass + FT(1)) / ρ / q_rai)^FT(1/(p_mass + 1))
end

"""
    terminal_velocity(q_rai, ρ)

where:
  - `q_rai` - rain water specific humidity
  - `ρ`     - density of air

Returns the mass weighted average rain terminal velocity assuming
Marshall Palmer 1948 distribution of rain drops.
"""
function terminal_velocity(q_rai::FT, ρ::FT) where {FT <: Real}

    return c_velo(ρ) * lambda(q_rai, ρ)^(-p_velo)
           * gamma(p_velo + p_mass + FT(1)) / gamma(p_mass + FT(1))

    # gamma(9/2)
    gamma_9_2 = FT(11.631728396567448)

    lambda::FT = (FT(8) * π * ρ_cloud_liq * MP_n_0 / ρ / q_rai)^FT(1 / 4)

    return gamma_9_2 * v_c / FT(6) * sqrt(grav / lambda)
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
    conv_q_liq_to_q_rai_accr(q_liq, q_rai, ρ)

where:
- `q_liq` - is the liquid water specific humidity
- `q_rai` - is the rain water specific humidity
- `ρ` - is the density of air

Returns the q_rai tendency due to collisions between cloud droplets
and rain drops (accretion) parametrized following Kessler 1995.
"""
function conv_q_liq_to_q_rai_accr(q_liq::FT, q_rai::FT, ρ::FT) where {FT<:Real}

    return MP_n_0 * c_area * c_velo * q_liq * E_col
           * gamma(p_area + p_velo + FT(1))
           * lambda(q_rai, ρ)^(-p_area - p_velo -FT(1))
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
function conv_q_rai_to_q_vap(qr::FT, q::PhasePartition{FT},
                             T::FT, p::FT, ρ::FT) where {FT<:Real}

    qv_sat = q_vap_saturation(T, ρ, q)
    q_v = q.tot - q.liq - q.ice
    S = q_v/qv_sat - FT(1)

    L = latent_heat_vapor(T)
    p_vs = saturation_vapor_pressure(T, Liquid())
    G::FT = FT(1) / (
              L / K_therm / T * (L / R_v / T - FT(1))
              + R_v * T / D_vapor / p_vs
            )

    return c_mass * p_mass * S * G * MP_n_0 / ρ / ρ_cloud_liq
           * (
               a_vent * lambda(q_rai, ρ)^(-p_mass + FT(1))
               * gamma(p_mass - FT(1))
               +
               b_vent * N_Sc^FT(1/3) * (FT(2) * c_velo / ν_air
               * lambda(q_rai, ρ)^FT(-2*p_mass - p_velo + 1))^FT(1/2)
               * gamma(FT(2*p_mass + p_velo - 1)/FT(2))
             )
end

end #module Microphysics.jl
