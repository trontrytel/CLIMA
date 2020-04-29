# Microphysics Module

The `Microphysics.jl` module describes a 1-moment bulk parameterization of cloud microphysical processes.
The module is based on the ideas of
[Kessler\_1995](https://www.sciencedirect.com/science/article/pii/016980959400090Z),
[Grabowski\_1998](https://journals.ametsoc.org/doi/full/10.1175/1520-0469%281998%29055%3C3283%3ATCRMOL%3E2.0.CO%3B2 )
and [Kaul\_et\_al\_2015](https://journals.ametsoc.org/doi/10.1175/MWR-D-14-00319.1)

The cloud microphysics variables are expressed as specific humidities:
  - `q_tot` - total water specific humidity,
  - `q_vap` - water vapor specific humidity,
  - `q_liq` - cloud water specific humidity,
  - `q_ice` - cloud ice specific humidity,
  - `q_rai` - rain specific humidity,
  - `q_sno` - snow specific humidity.

## Assumed particle size relationships

Particles are assumed to follow the mass(radius), cross section(radius) and terminal velocity(radius)
relationships defined as power laws.
The coefficients are defined in `MicrophysicsParameters` module and are shown in the Table below.
For rain and ice they correspond to spherical liquid water drops and ice particles, respectively.
There is no assumed particle shape for snow and the relationships are based on empirical values
from Grabowski 1998.

```math
m(r) = \alpha r^{\beta}
```
```math
a(r) = \gamma r^{\delta}
```
```math
v_{term}(r) = \zeta r^{\eta}
```
where:
 - ``r`` is the particle radius,
 - ``\alpha, \, \beta, \, \gamma, \, \delta, \, \zeta, \, \eta \,`` are the coefficients.

|    symbol          |         definition                              | units                    | default value                         |
|--------------------|-------------------------------------------------|--------------------------|---------------------------------------|
|``\alpha^{rai}``    | coefficient in mass(radius) for rain            | ``\frac{kg}{m^3}``       | ``\frac{4}{3} \, \pi \, \rho_{water}``|
|``\beta^{rai}``     | exponent in mass(radius) for rain               | -                        | ``3``                                 |
|``\gamma^{rai}``    | coefficient in cross section(radius) for rain   | -                        | ``\pi``                               |
|``\delta^{rai}``    | exponent in cross section(radius) for rain      | -                        | ``2``                                 |
|``\eta^{rai}``      | exponent in velocity(radius) for rain           | -                        | ``0.5``                               |
|                    |                                                 |                          |                                       |
|``\alpha^{ice}``    | coefficient in mass(radius) for ice             | ``\frac{kg}{m^3}``       | ``\frac{4}{3} \, \pi \, \rho_{ice}``  |
|``\beta^{ice}``     | exponent in mass(radius) for ice                | -                        | ``3``                                 |
|                    |                                                 |                          |                                       |
|``\alpha^{sno}``    | coefficient in mass(radius) for snow            | ``\frac{kg}{m^2}``       | ``0.1``                               |
|``\beta^{sno}``     | exponent in mass(radius) for snow               | -                        | ``2``                                 |
|``\gamma^{sno}``    | coefficient in cross section(radius) for snow   | -                        | ``0.3 \, \pi``                        |
|``\delta^{sno}``    | exponent in cross section(radius) for snow      | -                        | ``2``                                 |
|``\zeta^{sno}``     | coefficient in mass(radius) for snow            | ``\frac{m^{3/4}}{s}``    | ``2^{9/4}``                           |
|``\eta^{sno}``      | exponent in velocity(radius) for snow           | -                        | ``0.25``                              |

where:
 - ``\rho_{water}`` is the density of water,
 - ``\rho_{ice}`` is the density of ice.

The terminal velocity of an individual rain drop is defined by the balance between the gravitational acceleration
(taking into account the density difference between water and air) and the drag force.
Therefore the ``\zeta^{rai}`` is defined as
```math
\begin{equation}
\zeta^{rai} = \left(\frac{8}{3 \, C_{drag}} \left( \frac{\rho_{water}}{\rho} -1 \right) \right)^{1/2} (g)^{1/2}
\label{eq:vdrop}
\end{equation}
```
where:
 - ``g`` is the gravitational acceleration,
 - ``C_{drag}`` is the drag coefficient,
 - ``\rho`` is the density of air.\

---
**TODO**

It would be great to replace the above simple power laws with more accurate relationships.
For example:
[Khvorostyanov\_and\_Curry\_2002](https://journals.ametsoc.org/doi/pdf/10.1175/1520-0469%282002%29059%3C1872%3ATVODAC%3E2.0.CO%3B2)

---


## Assumed particle size distributions

The particle size distributions are assumed to follow Marshall-Palmer distribution (Marshall Palmer 1948 eq. 1):
```math
\begin{equation}
n(r) = n_{0} exp\left(- \lambda \, r \right)
\end{equation}
```
where:
 - ``n_{0}`` and ``\lambda`` are the Marshall-Palmer distribution parameters.

The ``n_0`` for rain and ice is assumed constant.
The ``n_0`` for snow is defined as
```math
\begin{equation}
n_0^{sno} = \mu^{sno} (\rho q_{sno})^{\nu^{sno}}
\end{equation}
```
where:
 - ``\mu^{sno}`` and ``\nu^{sno}`` are the coefficients

The coefficients are defined in `MicrophysicsParameters` module and are shown in the Table below.

|    symbol       |         definition                             | units              | default value       |
|-----------------|------------------------------------------------|--------------------|---------------------|
|``n_{0}^{rai}``  | rain drop size distribution parameter          | ``\frac{1}{m^4}``  | ``16 \cdot 10^6``   |
|``n_{0}^{ice}``  | cloud ice size distribution parameter          | ``\frac{1}{m^4}``  | ``2 \cdot 10^7``    |
|``\mu^{sno}``    | snow size distribution parameter coefficient   | ``?``              | ``4.36 \cdot 10^9`` |
|``\nu^{sno}``    | snow size distribution parameter exponent      | ``?``              | ``0.63``            |

The ``\lambda`` parameter is defined as
```math
\begin{equation}
\lambda = \left(\frac{\alpha \, n_0 \, \Gamma(\beta + 1)}{\rho \, q}\right)^{\frac{1}{\beta + 1}}
\end{equation}
```
where:
 - ``q`` is rain, snow or ice specific humidity
 - ``\alpha``, ``\beta`` and ``n_0`` are the corresponding mass-radius and size distribution parameters\
 - ``\Gamma()`` is the gamma function

The cloud-ice size distribution is only used when computing snow autoconversion rate.
In the derivation of different accretion rates the cloud ice, similar to cloud water,
is treated as continuous.

---
**TODO**

Do we want to keep the ``n_0`` for rain constant and ``n_0`` for snow empirical?

---

## Parameterized processes

Parameterized processes include:
  - diffusion of water vapour on cloud droplets and cloud ice crystals
    modeled as a relaxation to equilibrium
  - autoconversion,
  - accretion,
  - evaporation of rain water
  - melting - TODO
  - riming - TODO
  - sedimentation of rain and snow with mass weighted average terminal velocity
    (Water vapour, cloud water and cloud ice are part of the working fluid and do not sediment.)

Parameters used in the parameterization are defined in `MicrophysicsParameters` module. They consist of:

|    symbol                  |         definition                                        | units                    | default value         |
|----------------------------|-----------------------------------------------------------|--------------------------|-----------------------|
|``C_{drag}``                | rain drop drag coefficient                                | -                        | ``0.55``              |
|``\tau_{cond\_evap}``       | cloud water condensation/evaporation timescale            | ``s``                    | ``10``                |
|``\tau_{sub\_resub}``       | cloud ice condensation/evaporation timescale              | ``s``                    | ``10``                |
|``\tau_{acnv}``             | cloud to rain water autoconversion timescale              | ``s``                    | ``10^3``              |
|``q_{liq\_threshold}``      | cloud to rain water autoconversion threshold              | -                        | ``5 \cdot 10^{-4}``   |
|``r_{is}``                  | threshold particle radius between ice and snow            | ``m``                    | ``62.5 \cdot 10^{-6}``|
|``E_{cr}``                  | collision efficiency between rain drops and cloud droplets| -                        | ``0.8``               |
|``E_{cs}``                  | collision efficiency between snow and cloud droplets      | -                        | ``TODO``               |
|``E_{ir}``                  | collision efficiency between rain drops and cloud ice     | -                        | ``TODO``               |
|``E_{is}``                  | collision efficiency between snow and cloud ice           | -                        | ``TODO``               |
|``E_{rs}``                  | collision efficiency between rain drops and snow          | -                        | ``TODO``               |
|``a_{vent}, b_{vent}``      | rain drop ventilation factor coefficients                 | -                        | ``1.5 \;``,``\; 0.53``|
|``K_{therm}``               | thermal conductivity of air                               | ``\frac{J}{m \; s \; K}``| ``2.4 \cdot 10^{-2}`` |
|``\nu_{air}``               | kinematic viscosity of air                                | ``\frac{m^2}{s}``        | ``1.6 \cdot 10^{-5}`` |
|``D_{vapor}``               | diffusivity of water vapor                                | ``\frac{m^2}{s}``        | ``2.26 \cdot 10^{-5}``|

### Terminal velocity

The mass weighted terminal velocity ``v_t`` is defined following Ogura and Takahashi 1971
```math
\begin{equation}
v_t = \frac{\int_0^\infty n(r) \, m(r) \, v_{term}(r) \, dr}{\int_0^\infty n(r) \, m(r) \, dr}
\label{eq:vt}
\end{equation}
```
Integrating over the assumed Marshall-Palmer distribution results in and using the mass-radius and velocity-radius relationships
results in
```math
\begin{equation}
v_t = \zeta \lambda^{-\eta} \frac{\Gamma(\eta + \beta + 1)}{\Gamma(\beta + 1)}
\end{equation}
```
The default value of ``C_{drag}`` for rain is chosen such that the resulting mass averaged ``v_t`` is close to
the empirical terminal velocity formulation in Smolarkiewicz and Grabowski 1996.

---
**TODO**

Assuming a constant drag coefficient is an approximation and it should be size and flow dependent,
see [drag_coefficient](https://www.grc.nasa.gov/www/K-12/airplane/dragsphere.html).
In general it would be better to implement this:
[Khvorostyanov\_and\_Curry\_2002](https://journals.ametsoc.org/doi/pdf/10.1175/1520-0469%282002%29059%3C1872%3ATVODAC%3E2.0.CO%3B2)

---

### Cloud water condensation/evaporation

Condensation and evaporation of cloud water is parameterized as a relaxation to equilibrium value at the current time step.
```math
\begin{equation}
  \left. \frac{d \, q_{liq}}{dt} \right|_{cond, evap} = \frac{q^{eq}_{liq} - q_{liq}}{\tau_{cond\_evap}}
\end{equation}
```
where:
 - ``q^{eq}_{liq}`` - liquid water specific humidity in equilibrium,
 - ``q_{liq}`` - liquid water specific humidity,
 - ``\tau_{cond\_evap}`` - relaxation timescale (parameter in `MicrophysicsParameters` module).

### Cloud ice sublimation/resublimation

Sublimation and resublimation of cloud ice is parameterized as a relaxation to equilibrium value at the current time step.
```math
\begin{equation}
  \left. \frac{d \, q_{ice}}{dt} \right|_{sub, resub} = \frac{q^{eq}_{ice} - q_{ice}}{\tau_{sub\_resub}}
\end{equation}
```
where:
 - ``q^{eq}_{ice}`` - ice specific humidity in equilibrium,
 - ``q_{ice}`` - ice specific humidity,
 - ``\tau_{sub\_resub}`` - relaxation timescale (parameter in `MicrophysicsParameters` module).


---
**TODO**

Both ``\tau_{cond\_evap}`` and ``\tau_{sub\_resub}`` are assumed constant here.
It would be great to make the relaxation time a function of available condensation nuclei, turbulence intensity, etc.
See works by [prof Raymond Shaw](https://www.mtu.edu/physics/department/faculty/shaw/) for hints.

---

### Rain autoconversion

Rain autoconversion defines the rate of conversion form cloud to rain water due to collisions between cloud droplets. It is parameterized following Kessler 1995:
```math
\begin{equation}
  \left. \frac{d \, q_{rai}}{dt} \right|_{acnv} = \frac{max(0, q_{liq} - q_{liq\_threshold})}{\tau_{acnv}}
\end{equation}
```
where:
 - ``q_{liq}`` - liquid water specific humidity,
 - ``\tau_{acnv}`` - timescale (parameter in `MicrophysicsParameters` module),
 - ``q_{liq\_threshold}`` - autoconversion (parameter in `MicrophysicsParameters` module).

The default values of ``\tau_{acnv}`` and ``q_{liq\_threshold}`` are based on Smolarkiewicz and Grabowski 1996.

---
**TODO**

This is the simplest possible autoconversion parameterization.
It would be great to implement others and test the impact on precipitation.
See for example [Wood\_2005](https://journals.ametsoc.org/doi/full/10.1175/JAS3530.1)
Table 1 for other simple choices.

---

### Snow autoconversion

Snow autoconversion defines the rate of conversion form cloud ice to snow due
the growth of cloud ice by sublimation (deposition of water vapor).
It is defined as the change of mass of cloud ice for cloud ice particles
larger than threshold ``r_{is}``.

```math
\begin{equation}
  \left. \frac{d \, q_{sno}}{dt} \right|_{acnv} =
  \frac{1}{\rho} \frac{d}{dt} \left( \int_{r_{is}}^{\infty} m(r) n(r) dr \right) =
  \left. \frac{1}{\rho} \frac{dr}{dt} \right|_{r=r_{is}} m(r_{is}) n(r_{is}) + \frac{1}{\rho} \int_{r_{is}}^{\infty} \frac{dm}{dt} n(r) dr
\end{equation}
```
The ``\frac{dm}{dt}`` is obtained by solving the water vapor diffusion equation in spherical coordinates
and linking the changes in temperature at the drop surface to the changes in saturated vapor pressure via the Clausius-
Clapeyron equation, following Mason 1971.

For the simplest case of spherical particles and not taking into account ventilation effects:
```math
\begin{equation}
\frac{dm}{dt} = 4 \pi \, r \, (S - 1) \, G(T, p)
\label{eq:mass_rate}
\end{equation}
```
where:
 - ``S(q_{vap}, q_{vap}^{sat}) = \frac{q_{vap}}{q_{vap}^{sat}} - 1`` is commonly labeled as supersaturation,
 - ``q_{vap}^{sat}`` is the saturation vapor specific humidity,
 - ``G(T) = \left(\frac{L_s}{KT} \left(\frac{L_s}{R_v T} - 1 \right) + \frac{R_v T}{p_{vap}^{sat} D} \right)^{-1}``
     combines the effects of thermal conductivity and water diffusivity.
 - ``L_s`` is the latent heat of sublimation,
 - ``K_{thermo}`` is the thermal conductivity of air,
 - ``R_v`` is the gas constant of water vapor,
 - ``D_{vapor}`` is the diffusivity of water vapor

Using eq. (\ref{eq:mass_rate}) and the assumed mass-radius relationship we obtain
```math
\begin{equation}
\frac{dr}{dt} = \frac{4 \pi \, (S - 1) \, G(T, p)}{\alpha^{ice} \beta^{ice}} r^{2-\beta^{ice}}
\label{eq:r_rate}
\end{equation}
```
Finally the snow autoconversion rate is computed as
```math
\begin{equation}
  \left. \frac{d \, q_{sno}}{dt} \right|_{acnv} =
   \frac{1}{\rho} 4 \pi \, (S-1) \, G(T,p) \, n_0^{ice} \, exp(-\lambda_{ice} r_{is})
   \left( \frac{r_{is}^2}{\beta^{ice}} + \frac{r_{is} \lambda_{ice} +1}{\lambda_{ice}^2} \right)
\end{equation}
```

---
**TODO**

We should include ventialtion effects.

For non-spherical particles the mass rate of growth
should be multiplied by a function depending on the particle aspect ratio.
For functions proposed for different crystal habitats see
[Harrington_1995](https://journals.ametsoc.org/doi/pdf/10.1175/1520-0469%281995%29052%3C4344%3APOICCP%3E2.0.CO%3B2) Appentix B.

---

### Accretion

Accretion defines the rates of conversion between different categories due to collisions between particles.

For the case of collisions bewteen suspended water (cloud water or cloud ice)
and falling water (rain or snow) the sink of suspended water is defined as:

```math
\begin{equation}
\left. \frac{d \, q_{s}}{dt} \right|_{accr} =  - \int_0^\infty n_f(r) \, a_f(r) \, v_{term}(r) E_{sf} q_{s} dr
\label{eq:accr_1}
\end{equation}
```
where:
 - ``s`` subscript indicates suspended water category (cloud water or cloud ice)
 - ``f`` subscript indicates falling water category (rain or snow)
 - ``E_{sf}`` is the collision efficiency.

Integrating over the distribution yeilds:
```math
\begin{equation}
\left. \frac{d \, q_s}{dt} \right|_{accr} =
  - n_{0f} \, \zeta_f \, \gamma_f \, q_s \, E_{sf} \, \frac{\Gamma(\eta_f + \delta_f + 1)}{\lambda_f^{\eta_f + \delta_f +1}}
\label{eq:accr_fin}
\end{equation}
```
For the case of cloud water - rain, cloud ice  - snow collisions the sink of suspended water
becomes simply the source for falling water.
For the case of cloud water and snow collisions for temperatures below freezing
the sink of cloud water is a source for snow.
For temperatures above freezing the accreted cloud droplets along with some melted snow
are converted to rain.
In this case eq.(\ref{eq:accr_fin}) describes the sink of cloud water.
The sink of snow is proportional to the sink of cloud water with the coefficient ``\frac{c_w}{L_f}(T - T_{freeze})``.

The collisions between cloud ice - rain create snow.
The source of snow in this case is a sum of sinks from cloud ice and rain.
The sink of cloud ice is defined by eq. (\ref{eq:accr_1}).
The sink of rain is defined as:

```math
\begin{equation}
\left. \frac{d \, q_{rai}}{dt} \right|_{accr\_ri} =
  - \int_0^\infty \int_0^\infty
  \frac{1}{\rho} \, E_{ir} \, n_i(r_i) \, n_r(r_r) \, a_r(r_r) \, m_r(r_r)  \, v_{term}(r_r) \, d r_i d r_r
\label{eq:accr_ir}
\end{equation}
```
where:
 - ``E_{ir}`` is the collision efficiency between rain anc cloud ice
 - ``n_i`` and ``n_r`` are the cloud ice and rain size distributions
 - ``m_r``, ``a_r`` and ``v_{term}`` are the mass(radius), cross section(radius) and terminal velocity(radius) relations for rain
 - ``r_i`` and ``r_r`` mark integration over cloud ice and rain size distributions

Integrating eq.(\ref{eq:accr_ir}) yelds:
```math
\begin{equation}
\left. \frac{d \, q_{rai}}{dt} \right|_{accr\_ri} =
  - E_{ir} \, n_0^{rai} \, n_0^{ice} \, \alpha^{rai} \, \gamma^{rai} \, \zeta^{rai} \frac{1}{\lambda_{ice}}
  \frac{\Gamma(\beta^{rai} + \delta^{rai} + \eta^{rai} +1)}{\lambda_{rai}^{(\beta^{rai} + \delta^{rai} + \eta^{rai} + 1)}}
\end{equation}
```

Collisions between rain and snow result in snow in temperatures below freezing andin rain in temperatures above freezing.
The source term is defined as:
```math
\begin{equation}
\left. \frac{d \, q_i}{dt} \right|_{accr} =
    \int_0^\infty \int_0^\infty
    n_i(r_i) \, n_j(r_j) \, a(r_i, r_j) \, m_j(r_j) \, E_{ij} \, \left|v_{term}(r_i) - v_{term}(r_j)\right| \, dr_i dr_j
\label{eq:accr_sr1}
\end{equation}
```
where
- ``i`` stands for rain (``T>T_{freezing}``) or snow (``T<T_{freezing}``)
- ``j`` stands for snow (``T>T_{freezing}``) or rain (``T<T_{freezing}``)
- ``a(r_i, r_j)`` is the crossection of the two colliding particles

There are two additional assuptions that we make to integrate this:
- ``\left|v_{term}(r_i) - v_{term}(r_j)\right| \approx \left| v_{t_i} - v_{t_j} \right|``
  We approximate the teminal velocity difference for each particle pair with
  a difference between mass-weighted mean terminal velocityes and move it outside of the integral.
  See the discussion in [Ikawa\_and\_Saito\_1991](https://www.mri-jma.go.jp/Publish/Technical/DATA/VOL_28/28_005.pdf) page 88.

-  We assume that ``a(r_i, r_j) = \pi (r_i + r_j)^2``.
   This corresponds to a geometric formulation of the collision kernel, aka cylindrical fomulation,
   see [Wang\_et\_al\_2005](https://journals.ametsoc.org/doi/pdf/10.1175/JAS3655.1) for discussion.


The eq.(\ref{eq:accr_sr1}) can then be integrated as:
```math
\begin{align}
\left. \frac{d \, q_i}{dt} \right|_{accr} & =
    \pi \, n_0^{rai} \, n_0^{sno} \, \alpha_j \, E_{ij} \left| v_{t_i} - v_{t_j} \right|
    \int_0^\infty \int_0^\infty
    (r_i + r_j)^2
    r_{j}^{\beta_j} \,
    exp(- \lambda_j r_j) \,
    exp(- \lambda_i r_i) \,
    dr_i dr_j \\
    & =
    \pi \, n_0^{rai} \, n_0^{sno} \, \alpha_j \, E_{ij} \left| v_{t_i} - v_{t_j} \right|
    \left(
        \frac{2 \Gamma(\beta_j + 1)}{\lambda_s^3 \lambda_r^{\beta_r+1}}
        + \frac{2 \Gamma(\beta_j + 2)}{ \lambda_i^2 \lambda_j^{\beta_j+2}}
        + \frac{\Gamma(\beta_j + 3)}{\lambda_s \lambda_r^{\beta_r + 3}}
    \right)
\end{align}
```
---
**TODO**

Both of the assumptions needed to integrate the snow-rain accretion rate could be revisited:

The discussion on page 88 in
[Ikawa\_and\_Saito\_1991](https://www.mri-jma.go.jp/Publish/Technical/DATA/VOL_28/28_005.pdf)
suggests an alternative apporximation of the velocity difference.

The ``(r_i + r_j)^2`` assumption for the crossection is inconsistent with the
snow crossection used when modelling collisions with cloud water and cloud ice.

---


### Rain evaporation and snow sublimation

We start from a similar equation as when computing snow autoconversion rate
but integrate it from ``0`` to ``\infty``.

```math
\begin{equation}
  \left. \frac{dq}{dt} \right|_{evap\_subl} = \frac{1}{\rho} \int_{0}^{\infty} \frac{dm(r)}{dt} n(r) dr
\end{equation}
```
In contrast to eq.(\ref{eq:mass_rate}), now we are taking into account ventilation effects:

```math
\begin{equation}
  \frac{dm}{dt} = 4 \pi \, r \, (S - 1) \, G(T) \, F(r)
\label{eq:mass_rate2}
\end{equation}
```
where:
 - ``F(r)`` is the rain drop ventilation factor
 - supersaturation S is computed over water or ice
 - ``L`` is the latent heat of vaporization or sublimation.

Following Seifert and Beheng 2006 eq. 24 the ventilation factor is defined as:
```math
\begin{equation}
F(r) = a_{vent} + b_{vent}  N_{Sc}^{1/3} N_{Re}(r)^{1/2}
\label{eq:ventil_factor}
\end{equation}
```
where:
 - ``a_{vent}``, ``b_{vent}`` are coefficients,
 - ``N_{Sc}`` is the Schmidt number,
 - ``N_{Re}`` is the Reynolds number of a falling drop.
The Schmidt number is assumed constant:

```math
N_{Sc} = \frac{\nu_{air}}{D_{vapor}}
```
where:
 - ``\nu_{air}`` is the kinematic viscosity of air.
The Reynolds number of a spherical drop is defined as:

```math
N_{Re} = \frac{2 \, r \, v_{term}(r)}{\nu_{air}}
```
The final integral is:

```math
\begin{align}
\left. \frac{dq}{dt} \right|_{evap\_subl} & =
    \frac{4 \pi n_0}{\rho} (S - 1) G(T)
    \int_0^\infty
    \left(
       a_{vent} \, r +
       b_{vent} \left( \frac{\nu_{air}}{D_{vapor}} \right)^{\frac{1}{3}}
         \left( \frac{2 \zeta}{\nu_{air}} \right)^{\frac{1}{2}} r^{\frac{\eta + 3}{2}}
    \right)
    exp(-\lambda r) dr \\
    & =
    \frac{4 \pi n_0}{\rho} (S - 1) G(T) \lambda^{-2}
    \left(
       a_{vent} +
       b_{vent} \left( \frac{\nu_{air}}{D_{vapor}} \right)^{\frac{1}{3}}
         \left( \frac{2 \zeta}{\nu_{air}} \right)^{\frac{1}{2}}
         \Gamma \left( \frac{\eta + 5}{2} \right) \lambda^{- \frac{\eta +1}{2}}
    \right)
\end{align}
```

The values of ``a_{vent}`` and ``b_{vent}`` are chosen so that at ``q_{tot} = 15 g/kg`` and ``T=288K`` the resulting rain evaporation rate is close to the empirical rain evaporation rate from Smolarkiewicz and Grabowski 1996.

TODO - ``a_{vent}`` and ``b_{vent}`` for snow.
TODO - check units (with or without density?)

---
**TODO**

How to take into account the non-spherical snow shape? Modify the Reynolds number and growth equation.

---

### Snow melt

TODO - move ventilation factor definitions up top

If snow occurs in temperatures above freezing it is melting into rain.
The sink for rain is parameterized again as

```math
\begin{equation}
  \left. \frac{dq}{dt} \right|_{melt} = \frac{1}{\rho} \int_{0}^{\infty} \frac{dm(r)}{dt} n(r) dr
\end{equation}
```
For snow melt
```math
\begin{equation}
  \frac{dm}{dt} = 4 \pi \, r \, \frac{K_{thermo}}{L_f} (T - T_{freeze}) \, F(r)
\label{eq:mass_rate3}
\end{equation}
```
where:
 - ``F(r)`` is the ventilation factor defined as in \ref{eq:ventil_factor}
 - ``L_f`` is the latent heat of freezing.

```math
\begin{equation}
\left. \frac{dq}{dt} \right|_{evap\_subl} =
    \frac{4 \pi \, n_0 \, K_{thermo}}{\rho \, L_f} (T - T_{freeze}) \lambda^{-2}
    \left(
       a_{vent} +
       b_{vent} \left( \frac{\nu_{air}}{D_{vapor}} \right)^{\frac{1}{3}}
         \left( \frac{2 \zeta}{\nu_{air}} \right)^{\frac{1}{2}}
         \Gamma \left( \frac{\eta + 5}{2} \right) \lambda^{- \frac{\eta +1}{2}}
    \right)
\end{equation}
```

## Examples (TODO)

```@example rain_terminal_velocity
using CLIMA.Microphysics
using Plots

# eq. 5d in Smolarkiewicz and Grabowski 1996
# https://doi.org/10.1175/1520-0493(1996)124<0487:TTLSLM>2.0.CO;2
function terminal_velocity_empirical(q_rai::DT, q_tot::DT, ρ::DT, ρ_air_ground::DT) where {DT<:Real}
    rr  = q_rai / (DT(1) - q_tot)
    vel = DT(14.34) * ρ_air_ground^DT(0.5) * ρ^-DT(0.3654) * rr^DT(0.1346)
    return vel
end

q_rain_range = range(1e-8, stop=5e-3, length=100)
ρ_air, q_tot, ρ_air_ground = 1.2, 20 * 1e-3, 1.22

plot(q_rain_range * 1e3,  [terminal_velocity(q_rai, ρ_air) for q_rai in q_rain_range], xlabel="q_rain [g/kg]", ylabel="velocity [m/s]", title="Average terminal velocity of rain", label="CLIMA")
plot!(q_rain_range * 1e3, [terminal_velocity_empirical(q_rai, q_tot, ρ_air, ρ_air_ground) for q_rai in q_rain_range], label="Empirical")
savefig("rain_terminal_velocity.svg") # hide
nothing # hide
```
![](rain_terminal_velocity.svg)


The default value of collision efficiency ``E_{coll}`` is set to 0.8 so that
the resulting accretion rate is close to the empirical accretion rate in Smolarkiewicz and Grabowski 1996.
Assuming a constant ``E_{col}`` is an approximation,
see for example [collision efficiency](https://journals.ametsoc.org/doi/10.1175/1520-0469%282001%29058%3C0742%3ACEODIA%3E2.0.CO%3B2).

```@example accretion
using CLIMA.Microphysics
using Plots

# eq. 5b in Smolarkiewicz and Grabowski 1996
# https://doi.org/10.1175/1520-0493(1996)124<0487:TTLSLM>2.0.CO;2
function accretion_empirical(q_rai::DT, q_liq::DT, q_tot::DT) where {DT<:Real}
    rr  = q_rai / (DT(1) - q_tot)
    rl  = q_liq / (DT(1) - q_tot)
    return DT(2.2) * rl * rr^DT(7/8)
end

# some example values
q_rain_range = range(1e-8, stop=5e-3, length=100)
ρ_air, q_liq, q_tot = 1.2, 5e-4, 20e-3

plot(q_rain_range * 1e3,  [conv_q_liq_to_q_rai_accr(q_liq, q_rai, ρ_air) for q_rai in q_rain_range], xlabel="q_rain [g/kg]", ylabel="accretion rate [1/s]", title="Accretion", label="CLIMA")
plot!(q_rain_range * 1e3, [accretion_empirical(q_rai, q_liq, q_tot) for q_rai in q_rain_range], label="empirical")
savefig("accretion_rate.svg") # hide
nothing # hide
```
![](accretion_rate.svg)



```@example rain_evaporation
using CLIMA.Microphysics
using CLIMA.MoistThermodynamics

using CLIMAParameters
using CLIMAParameters.Planet: R_d, planet_radius, grav, MSLP
struct EarthParameterSet <: AbstractEarthParameterSet end
const param_set = EarthParameterSet()

using Plots

# eq. 5c in Smolarkiewicz and Grabowski 1996
# https://doi.org/10.1175/1520-0493(1996)124<0487:TTLSLM>2.0.CO;2
function rain_evap_empirical(q_rai::DT, q::PhasePartition, T::DT, p::DT, ρ::DT) where {DT<:Real}

    q_sat  = q_vap_saturation(param_set, T, ρ, q)
    q_vap  = q.tot - q.liq
    rr     = q_rai / (DT(1) - q.tot)
    rv_sat = q_sat / (DT(1) - q.tot)
    S      = q_vap/q_sat - DT(1)

    ag, bg = 5.4 * 1e2, 2.55 * 1e5
    G = DT(1) / (ag + bg / p / rv_sat) / ρ

    av, bv = 1.6, 124.9
    F = av * (ρ/DT(1e3))^DT(0.525)  * rr^DT(0.525) + bv * (ρ/DT(1e3))^DT(0.7296) * rr^DT(0.7296)

    return DT(1) / (DT(1) - q.tot) * S * F * G
end

# example values
T, p = 273.15 + 15, 90000.
ϵ = 1. / molmass_ratio
p_sat = saturation_vapor_pressure(param_set, T, Liquid())
q_sat = ϵ * p_sat / (p + p_sat * (ϵ - 1.))
q_rain_range = range(1e-8, stop=5e-3, length=100)
q_tot = 15e-3
q_vap = 0.15 * q_sat
q_ice = 0.
q_liq = q_tot - q_vap - q_ice
q = PhasePartition(q_tot, q_liq, q_ice)
R = gas_constant_air(param_set, q)
ρ = p / R / T

plot(q_rain_range * 1e3,  [conv_q_rai_to_q_vap(q_rai, q, T, p, ρ) for q_rai in q_rain_range], xlabel="q_rain [g/kg]", ylabel="rain evaporation rate [1/s]", title="Rain evaporation", label="CLIMA")
plot!(q_rain_range * 1e3, [rain_evap_empirical(q_rai, q, T, p, ρ) for q_rai in q_rain_range], label="empirical")
savefig("rain_evaporation_rate.svg") # hide
nothing # hide
```
![](rain_evaporation_rate.svg)



```@meta
CurrentModule = CLIMA.Microphysics
```

## User interface

```@docs
terminal_velocity
conv_q_vap_to_q_liq
conv_q_liq_to_q_rai_acnv
conv_q_liq_to_q_rai_accr
conv_q_rai_to_q_vap
```

## References

@article{Grabowski_and_Smolarkiewicz_1996,
author = {Grabowski, Wojciech W. and Smolarkiewicz, Piotr K.},
title = {Two-Time-Level Semi-Lagrangian Modeling of Precipitating Clouds},
journal = {Monthly Weather Review},
volume = {124},
number = {3},
pages = {487-497},
year = {1996},
doi = {10.1175/1520-0493(1996)124<0487:TTLSLM>2.0.CO;2}}

@article{Grabowski_1998,
author = {Grabowski, Wojciech W.},
title = {Toward Cloud Resolving Modeling of Large-Scale Tropical Circulations: A Simple Cloud Microphysics Parameterization},
journal = {Journal of the Atmospheric Sciences},
volume = {55},
number = {21},
pages = {3283-3298},
year = {1998},
doi = {10.1175/1520-0469(1998)055<3283:TCRMOL>2.0.CO;2}}

@article{Kaul_et_al_2015,
author = {Kaul, Colleen M. and Teixeira, João and Suzuki, Kentaroh},
title = {Sensitivities in Large-Eddy Simulations of Mixed-Phase Arctic Stratocumulus Clouds Using a Simple Microphysics Approach},
journal = {Monthly Weather Review},
volume = {143},
number = {11},
pages = {4393-4421},
year = {2015},
doi = {10.1175/MWR-D-14-00319.1}}

@article{Kessler_1995,
author = {Kessler, E.},
title = {On the continuity and distribution of water substance in atmospheric circulations},
journal = {Atmospheric Research},
volume = {38},
number = {1},
pages = {109 - 145},
year = {1995},
doi = {10.1016/0169-8095(94)00090-Z}}

@book{Mason_1971,
author = {Mason, B. J.},
title = {The Physics of Clouds},
publisher = {Oxford Univ. Press},
year = {1971}}

@article{Marshall_and_Palmer_1948,
author = {Marshall, J. S. and Palmer, W. Mc K.},
title = {The distribution of raindrops with size},
journal = {Journal of Meteorology},
volume = {5},
number = {4},
pages = {165-166},
year = {1948},
doi = {10.1175/1520-0469(1948)005<0165:TDORWS>2.0.CO;2}}

@article{Ogura_and_Takahashi_1971,
author = {Oqura, Yoshimitsu and Takahashi, Tsutomu},
title = {Numerical simulation of the life cycle of a thunderstorm cell},
journal = {Monthly Weather Review},
volume = {99},
number = {12},
pages = {895-911},
year = {1971},
doi = {10.1175/1520-0493(1971)099<0895:NSOTLC>2.3.CO;2}}

@article{Seifert_and_Beheng_2006,
author={Seifert, A. and Beheng, K. D.},
title={A two-moment cloud microphysics parameterization for mixed-phase clouds. Part 1: Model description},
journal={Meteorology and Atmospheric Physics},
year={2006},
volume={92},
number={1},
pages={45--66},
doi={10.1007/s00703-005-0112-4}}
