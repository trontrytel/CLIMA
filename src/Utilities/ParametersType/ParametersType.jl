module ParametersType
using Unitful
import Unitful: AbstractQuantity, Units

export @parameter, @exportparameter, ParametersType

"""
    Parameter{sym} <: Base.AbstractIrrational

Number type representing a constant parameter value denoted by the symbol `sym`.

!!! note

AbstractIrrational is used here inorder to inherit the behavior from the Base.
ParametersType need not be Irrational numbers.

"""
struct Parameter{sym} <: Base.AbstractIrrational end

Base.show(io::IO, x::Parameter{S}) where {S} = print(io, "$(string(x))")

Base.:(==)(x::Parameter, y::AbstractQuantity) = (getval(x) == y)
Base.:(==)(x::Parameter, y::AbstractFloat) = (getval(x) == y)
Base.:(==)(x::Parameter, y::Irrational) = (getval(x) == y)
Base.:(==)(x::Parameter, y::Rational) = (getval(x) == y)
Base.:(==)(x::AbstractQuantity, y::Parameter) = (x == getval(y))
Base.:(==)(x::AbstractFloat, y::Parameter) = (x == getval(y))
Base.:(==)(x::Irrational, y::Parameter) = (x == getval(y))
Base.:(==)(x::Rational, y::Parameter) = (x == getval(y))
Base.:(==)(x::Parameter, y::Parameter) = (getval(x) == getval(y))
Base.:<(x::Parameter, y::Parameter) = (getval(x) < getval(y))
Base.:<=(x::Parameter, y::Parameter) = (getval(x) <= getval(y))
Base.:*(x::Parameter, y::Parameter) = getval(x) * getval(y)
Base.:*(x::Parameter, y::AbstractQuantity) = getval(x) * y
Base.:*(x::AbstractQuantity, y::Parameter) = x * getval(y)
Base.:*(x::Parameter, y::Units) = getval(x) * y
Base.:*(x::Units, y::Parameter) = x * getval(x)
Base.:/(x::Parameter, y::Parameter) = getval(x) / getval(y)
Base.:/(x::Parameter, y::AbstractQuantity) = getval(x) / y
Base.:/(x::AbstractQuantity, y::Parameter) = x / getval(y)
Base.:/(x::Parameter, y::Units) = getval(x) / y
Base.:/(x::Units, y::Parameter) = x / getval(y)
Base.:+(x::Parameter, y::Parameter) = getval(x) + getval(y)
Base.:+(x::Parameter, y::AbstractQuantity) = getval(x) + y
Base.:+(x::AbstractQuantity, y::Parameter) = x + getval(y)
Base.:+(x::Parameter, y::Units) = getval(x) + y
Base.:+(x::Units, y::Parameter) = x + getval(x)
Base.:-(x::Parameter, y::Parameter) = getval(x) - getval(y)
Base.:-(x::Parameter, y::AbstractQuantity) = getval(x) - y
Base.:-(x::AbstractQuantity, y::Parameter) = x - getval(y)
Base.:-(x::Parameter, y::Units) = getval(x) - y
Base.:-(x::Units, y::Parameter) = x - getval(y)
Base.hash(x::Parameter, h::UInt) = 3*objectid(x) - h
Base.widen(::Type{T}) where {T<:Parameter} = T
Base.round(x::Parameter, r::RoundingMode) = round(float(x), r)
getval() = nothing
getval(x::Number) = x

"""
    @parameter sym val desc doexport=false
    @parameter(sym, val, desc, doexport-false)

Define a new `Parameter` value, `sym`, with value `val` and description string
`desc`. If `doexport == true` then `sym` is exported, e.g., the command
`export \$sym` is added
"""
macro parameter(sym, val, desc, doexport=false)
  esym = esc(sym)
  qsym = esc(Expr(:quote, sym))
  ev = getval(@eval(__module__, $val))

  exportcmd = doexport ? :(export $sym) : ()

  quote
    $exportcmd
    const $esym = Parameter{$qsym}()
    Base.Float64(::Parameter{$qsym}) = $(Float64(ustrip(ev)) * upreferred(unit(ev)))
    Base.Float32(::Parameter{$qsym}) = $(Float32(ustrip(ev)) * upreferred(unit(ev)))
    Base.string(::Parameter{$qsym}) = $(string(ev))
    ParametersType.getval(::Parameter{$qsym}) = $(esc(ev))
    """
        $($qsym)

    $($desc)

    # Examples
    ```
    julia> $($qsym)
    $($(string(ev)))
    ```
    """
    $sym
  end
end

global symbols = []

# TODO: figure out how to get this to call @parameter with doexport=true
macro exportparameter(sym, val, desc)
  esym = esc(sym)
  qsym = esc(Expr(:quote, sym))
  ev = getval(@eval(__module__, $val))

  exportcmd = :(export $sym)

  quote
    $exportcmd
    const $esym = Parameter{$qsym}()
    if !@isdefined symbols
        global symbols = []
    end
    push!(symbols, $esym)
    # FIXME: shouldn't really allow constructors to return other types
    Base.Float64(::Parameter{$qsym}) = $(Float64(ustrip(ev)) * upreferred(unit(ev)))
    Base.Float32(::Parameter{$qsym}) = $(Float32(ustrip(ev)) * upreferred(unit(ev)))
    Base.string(::Parameter{$qsym}) = $(string(ev))
    ParametersType.getval(::Parameter{$qsym}) = $(esc(upreferred(getval(ev))))
    """
        $($qsym)

    $($desc)

    # Examples
    ```
    julia> $($qsym)
    $($(string(ev)))
    ```
    """
    $sym
  end
end

end
