"""
    MicrophysicsParameters

Module containing 1-moment bulk microphysics parameters.
"""
module MicrophysicsParameters
using CLIMA.ParametersType
using CLIMA.PlanetParameters

# rain: liquid and spherical
# ice:  ice and spherical
# snow:

# shape relationships
# mass(radius)     = c_mass      * radius^p_mass
# velocity(radius) = c_velo(rho) * radius^p_velo
# area(radius)     = c_area      * radius^p_area

# assumed size distribution
# n(r) = MP_n_0 * exp (-lambda r)

@exportparameter c_mass_rai 4/3. * π * ρ_cloud_liq "coefficient in mass(radius)"
@exportparameter p_mass_rai 3 "exponent in mass(radius)"
@exportparameter p_velo_rai 0.5 "exponent in velocity(radius)"
@exportparameter c_area_rai π "coefficient in area(radius) relation"
@exportparameter p_area_rai 2 "exponent in area(radius) relation"

@exportparameter c_mass_ice 4/3. * π * ρ_cloud_ice "coefficient in mass(radius)"
@exportparameter p_mass_ice 3 "exponent in mass(radius)"

@exportparameter c_mass_sno 1e-1 "coefficient in mass(radius)"
@exportparameter p_mass_sno 2 "exponent in mass(radius)"
@exportparameter c_velo_sno 2^(9/4) "coefficient in velocity(radius)"
@exportparameter p_velo_sno 0.25 "exponent in velocity(radius)"

@exportparameter MP_n_0_rai 8e6 * 2 "Marshall-Palmer distribution n_0 coeff for rain [1/m4]"
@exportparameter MP_n_0_ice 1e7 * 2 "Marshall-Palmer distribution n_0 coeff for ice [1/m4]"
@exportparameter MP_sno_c 4.36 * 1e9 "Marshall-Palmer distribution coefficient tocompute n_0 for snow [?]"
@exportparameter MP_sno_p 0.63 "Marshall-Palmer distribution coefficient to compute n_0 for snow [?]"

@exportparameter r_ice_snow 62.5 * 1e-6 "threshold between ice and snow"

@exportparameter τ_cond_evap_liq 10 "condensation/evaporation timescale [s]"
@exportparameter τ_cond_evap_ice 10 "sublimation/resublimation timescale [s]"

@exportparameter τ_acnv 1e3 "autoconversion timescale [s]  ∈(1e3, 1e4) "
@exportparameter q_liq_threshold 5e-4 "autoconv. threshold [-] ∈(.5, 1) * 1e-3"

@exportparameter E_col_rai 0.8 "collision efficiency for rain [-]"
@exportparameter E_col_sno 0.2 "collision efficiency for snow [-]" #TODO or 0.8?
@exportparameter α_geom_rai 1 "geometric factor for rain [-]"
@exportparameter α_geom_sno 0.3 "geometric factor for rain [-]"

@exportparameter C_drag 0.55 "drag coefficient for rain drops [-]"

@exportparameter a_vent 1.5 "ventilation factor coefficient [-]"
@exportparameter b_vent 0.53 "ventilation factor coefficient [-]"
@exportparameter K_therm 2.4e-2 "thermal conductivity of air [J/m/s/K] "
@exportparameter D_vapor 2.26e-5 "diffusivity of water vapor [m2/s]"
@exportparameter ν_air 1.6e-5 "kinematic viscosity of air [m2/s]"
@exportparameter N_Sc ν_air/D_vapor "Schmidt number [-]"
end
