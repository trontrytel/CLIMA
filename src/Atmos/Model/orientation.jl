# TODO: add Coriolis vectors
import ..PlanetParameters: grav, planet_radius
using ..DGmethods: space_unit, gravp_unit, acc_unit
export Orientation, NoOrientation, FlatOrientation, SphericalOrientation

abstract type Orientation
end


function vars_aux(m::Orientation, T)
  @vars begin
    Φ::U(T, u"J/kg") # gravitational potential
    ∇Φ::SVector{3, U(T, u"J/kg/m")}
  end
end

gravitational_potential(::Orientation, aux::Vars) = aux.orientation.Φ
∇gravitational_potential(::Orientation, aux::Vars) = aux.orientation.∇Φ
altitude(orientation::Orientation, aux::Vars) = gravitational_potential(orientation, aux) / grav
vertical_unit_vector(orientation::Orientation, aux::Vars) = ∇gravitational_potential(orientation, aux) / grav


"""
    NoOrientation

No gravitional force or potential.
"""
struct NoOrientation <: Orientation
end
function vars_aux(m::NoOrientation, T)
  @vars()
end
atmos_init_aux!(::NoOrientation, ::AtmosModel, aux::Vars, geom::LocalGeometry) = nothing
gravitational_potential(::NoOrientation, aux::Vars) = -zero(eltype(aux))
∇gravitational_potential(::NoOrientation, aux::Vars) = SVector{3,eltype(aux)}(0,0,0)
altitude(orientation::NoOrientation, aux::Vars) = -zero(eltype(aux))

"""
    SphericalOrientation <: Orientation

Gravity acts towards the origin `(0,0,0)`, and the gravitational potential is relative
to the surface of the planet.
"""
struct SphericalOrientation <: Orientation
end
function atmos_init_aux!(::SphericalOrientation, ::AtmosModel, aux::Vars, geom::LocalGeometry)
  normcoord = norm(aux.coord)
  aux.orientation.Φ = grav * (normcoord - planet_radius)
  aux.orientation.∇Φ = grav / normcoord .* aux.coord
end

"""
    FlatOrientation <: Orientation

Gravity acts in the third coordinate, and the gravitational potential is relative to
`coord[3] == 0`.
"""
struct FlatOrientation <: Orientation
  # for Coriolis we could add latitude?
end
function atmos_init_aux!(::FlatOrientation, m::AtmosModel, aux::Vars, geom::LocalGeometry)
  aux.orientation.Φ = grav * aux.coord[3]
  u_grad_pot = gravp_unit(m) / space_unit(m)
  v = (0, 0, grav / acc_unit(m)).*u_grad_pot
  aux.orientation.∇Φ = SVector{3, U(eltype(aux), u_grad_pot)}(v)
end
