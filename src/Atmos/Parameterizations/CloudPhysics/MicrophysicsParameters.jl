"""
    MicrophysicsParameters

Module containing 1-moment bulk microphysics parameters.
"""
module MicrophysicsParameters
using CLIMA.ParametersType
using CLIMA.PlanetParameters

# liquid spherical rain drops:
# mass(radius)     = c_mass      * r^p_mass
# velocity(radius) = c_velo(rho) * r^p_velo
# area(radius)     = c_area      * r^p_area

@exportparameter c_mass 4/3. * π * ρ_cloud_liq "coefficient in mass(radius)"
@exportparameter p_mass 3 "exponent in mass(radius)"
@exportparameter p_velo 0.5 "exponent in velocity(radius)"
@exportparameter c_area π "coefficient in area(radius) relation"
@exportparameter p_area 2 "exponent in area(radius) relation"

@exportparameter MP_n_0 8e6 * 2 "Marshall-Palmer distribution n_0 coeff [1/m4]"
@exportparameter C_drag 0.55 "drag coefficient for rain drops [-]"

@exportparameter τ_cond_evap_liq 10 "condensation/evaporation timescale [s]"
@exportparameter τ_cond_evap_ice 10 "sublimation/resublimation timescale [s]"

@exportparameter τ_acnv 1e3 "autoconversion timescale [s]  ∈(1e3, 1e4) "
@exportparameter q_liq_threshold 5e-4 "autoconv. threshold [-] ∈(.5, 1) * 1e-3"

@exportparameter E_col 0.8 "collision efficiency [-]"

@exportparameter a_vent 1.5 "ventilation factor coefficient [-]"
@exportparameter b_vent 0.53 "ventilation factor coefficient [-]"
@exportparameter K_therm 2.4e-2 "thermal conductivity of air [J/m/s/K] "
@exportparameter D_vapor 2.26e-5 "diffusivity of water vapor [m2/s]"
@exportparameter ν_air 1.6e-5 "kinematic viscosity of air [m2/s]"
@exportparameter N_Sc ν_air/D_vapor "Schmidt number [-]"

end
