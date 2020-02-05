#### Subsidence model
export AbstractSubsidence,
       NoSubsidence,
       ConstantSubsidence,
       subsidence_velocity

abstract type AbstractSubsidence{FT<:AbstractFloat} end

struct NoSubsidence{FT} <: AbstractSubsidence{FT} end

struct ConstantSubsidence{FT} <: AbstractSubsidence{FT}
  D::U(FT,:frequency)
end

subsidence_velocity(::NoSubsidence{FT}, z::U(FT,:space)) where {FT} = FT(0) * unit_alias(:velocity) #FIXME
subsidence_velocity(subsidence::ConstantSubsidence{FT}, z::U(FT,:space)) where {FT} = subsidence.D*z
