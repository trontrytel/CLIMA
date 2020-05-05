using Test
using CLIMA.Microphysics
using CLIMA.MicrophysicsParameters
using CLIMA.MoistThermodynamics
using CLIMA.PlanetParameters: molmass_ratio, R_v

using CLIMA.ParametersType #TODO - ugly! needed to unpack parameters

@testset "RainFallSpeed" begin

    # eq. 5d in Smolarkiewicz and Grabowski 1996
    # https://doi.org/10.1175/1520-0493(1996)124<0487:TTLSLM>2.0.CO;2
    function terminal_velocity_empir(
        q_rai::FT,
        q_tot::FT,
        ρ::FT,
        ρ_air_ground::FT,
    ) where {FT <: Real}
        rr = q_rai / (1 - q_tot)
        vel = FT(14.34) * ρ_air_ground^FT(0.5) * ρ^-FT(0.3654) * rr^FT(0.1346)
        return vel
    end

    # some example values
    q_rain_range = range(1e-8, stop = 5e-3, length = 10)
    ρ_air, q_tot, ρ_air_ground = 1.2, 20 * 1e-3, 1.22

    #TODO - ugly way to unpack parameters
    n0_rai_ = ParametersType.getval(n0_rai)
    α_rai_ = ParametersType.getval(α_rai)
    β_rai_ = ParametersType.getval(β_rai)
    η_rai_ = ParametersType.getval(η_rai)
    ζ_rai_ = ζ_rai(ρ_air)

    for q_rai in q_rain_range

        λ_rai_ = lambda(q_rai, ρ_air, n0_rai_, α_rai_, β_rai_)

        @test terminal_velocity(q_rai, λ_rai_, β_rai_, ζ_rai_, η_rai_) ≈
              terminal_velocity_empir(q_rai, q_tot, ρ_air, ρ_air_ground) atol =
            0.2 * terminal_velocity_empir(q_rai, q_tot, ρ_air, ρ_air_ground)

    end
end

@testset "CloudCondEvap" begin

    q_liq_sat = 5e-3
    frac = [0.0, 0.5, 1.0, 1.5]

    #TODO - ugly way to unpack parameters
    τ_cond_evap_ = ParametersType.getval(τ_cond_evap)

    for fr in frac
        q_liq = q_liq_sat * fr

        @test conv_q_vap_to_q_liq_ice(
                  PhasePartition(0.0, q_liq_sat, 0.0),
                  PhasePartition(0.0, q_liq, 0.0),
                  τ_cond_evap_,
                  Liquid(),
            ) ≈ (1 - fr) * q_liq_sat / τ_cond_evap_
    end
end

@testset "RainAutoconversion" begin

    #TODO - ugly way to unpack parameters
    τ_acnv_ = ParametersType.getval(τ_acnv)

    q_liq_small = 0.5 * q_liq_threshold
    @test conv_q_liq_to_q_rai(q_liq_small, τ_acnv_) == 0.0

    q_liq_big = 1.5 * q_liq_threshold
    @test conv_q_liq_to_q_rai(q_liq_big, τ_acnv_) == 0.5 * q_liq_threshold / τ_acnv_
end

@testset "RainAccretion" begin

    # eq. 5b in Smolarkiewicz and Grabowski 1996
    # https://doi.org/10.1175/1520-0493(1996)124<0487:TTLSLM>2.0.CO;2
    function accretion_empir(q_rai::FT, q_liq::FT, q_tot::FT) where {FT <: Real}
        rr = q_rai / (FT(1) - q_tot)
        rl = q_liq / (FT(1) - q_tot)
        return FT(2.2) * rl * rr^FT(7 / 8)
    end

    # some example values
    q_rain_range = range(1e-8, stop = 5e-3, length = 10)
    ρ_air, q_liq, q_tot = 1.2, 5e-4, 20e-3

    #TODO - ugly way to unpack parameters
    n0_rai_ = ParametersType.getval(n0_rai)
    α_rai_ = ParametersType.getval(α_rai)
    β_rai_ = ParametersType.getval(β_rai)
    γ_rai_ = ParametersType.getval(γ_rai)
    δ_rai_ = ParametersType.getval(δ_rai)
    ζ_rai_ = ζ_rai(ρ_air)
    η_rai_ = ParametersType.getval(η_rai)
    E_lr_ = ParametersType.getval(E_lr)

    for q_rai in q_rain_range

        λ_rai_ = lambda(q_rai, ρ_air, n0_rai_, α_rai_, β_rai_)

        @test accretion(q_liq, q_rai, E_lr_, n0_rai_, λ_rai_, γ_rai_,
                δ_rai_, ζ_rai_, η_rai_) ≈ accretion_empir(q_rai, q_liq,
                q_tot) atol = (0.1 * accretion_empir(q_rai, q_liq, q_tot))
    end
end

@testset "RainEvaporation" begin

    # eq. 5c in Smolarkiewicz and Grabowski 1996
    # https://doi.org/10.1175/1520-0493(1996)124<0487:TTLSLM>2.0.CO;2
    function rain_evap_empir(
        q_rai::FT,
        q::PhasePartition,
        T::FT,
        p::FT,
        ρ::FT,
    ) where {FT <: Real}

        q_sat = q_vap_saturation(T, ρ, q)
        q_vap = q.tot - q.liq
        rr = q_rai / (1 - q.tot)
        rv_sat = q_sat / (1 - q.tot)
        S = q_vap / q_sat - 1

        ag, bg = FT(5.4 * 1e2), FT(2.55 * 1e5)
        G = FT(1) / (ag + bg / p / rv_sat) / ρ

        av, bv = FT(1.6), FT(124.9)
        F =
            av * (ρ / FT(1e3))^FT(0.525) * rr^FT(0.525) +
            bv * (ρ / FT(1e3))^FT(0.7296) * rr^FT(0.7296)

        return 1 / (1 - q.tot) * S * F * G
    end

    # example values
    T, p = 273.15 + 15, 90000.0
    ϵ = 1.0 / molmass_ratio
    p_sat = saturation_vapor_pressure(T, Liquid())
    q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.0))
    q_rain_range = range(1e-8, stop = 5e-3, length = 10)
    q_tot = 15e-3
    q_vap = 0.15 * q_sat
    q_ice = 0.0
    q_liq = q_tot - q_vap - q_ice
    q = PhasePartition(q_tot, q_liq, q_ice)
    R = gas_constant_air(q)
    ρ = p / R / T

    #TODO - ugly way to unpack parameters
    n0_rai_ = ParametersType.getval(n0_rai)
    α_rai_ = ParametersType.getval(α_rai)
    β_rai_ = ParametersType.getval(β_rai)
    γ_rai_ = ParametersType.getval(γ_rai)
    δ_rai_ = ParametersType.getval(δ_rai)
    ζ_rai_ = ζ_rai(ρ)
    η_rai_ = ParametersType.getval(η_rai)
    E_lr_ = ParametersType.getval(E_lr)
    a_vent_rai_ = ParametersType.getval(a_vent_rai)
    b_vent_rai_ = ParametersType.getval(b_vent_rai)
    K_therm_ = ParametersType.getval(K_therm)
    R_v_ = ParametersType.getval(R_v)
    ν_air_ = ParametersType.getval(ν_air)
    D_vapor_ = ParametersType.getval(D_vapor)
    S_liq_ = supersaturation(q, ρ, T, Liquid())
    G_liq_ = G_func(T, K_therm_, R_v_, D_vapor_, Liquid())

    for q_rai in q_rain_range

        λ_rai_ = lambda(q_rai, ρ, n0_rai_, α_rai_, β_rai_)

        @test evaporation_sublimation(q_rai, ρ, S_liq_, G_liq_, a_vent_rai_,
                b_vent_rai_, ν_air_, D_vapor_, n0_rai_, λ_rai_, ζ_rai_,
                η_rai_) ≈ rain_evap_empir(q_rai, q, T, p, ρ) atol = -0.5 *
                rain_evap_empir(q_rai, q, T, p, ρ)
    end
end
