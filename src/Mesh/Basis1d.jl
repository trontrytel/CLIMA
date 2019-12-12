module Basis1d
import GaussQuadrature
using LinearAlgebra

export B1d, jacobip

struct B1d{T <:Integer, FT <:AbstractFloat}

    qm1::T               # Number of quadrature 1D points on mesh 1

    m1_r::Vector{FT}     # 1D GLL points for mesh 1
    m1_w::Vector{FT}     # 1D GLL weights for mesh 1 

    m1_phi::Array{FT,2}  # Legendre polynomials on mesh 1 

    trans::Array{FT,2}   # transfer from Lagrangian to  Legendre basis (u_Leg = trans * u_Lag)
#-----------------------------------------
function B1d(deg::T, ::Type{FT}) where {T <:Integer, FT <:AbstractFloat} # inner constructor

   alpha, beta = 0, 0 # Legendre polynomials
   qm1 = deg + T(1)
   m1_r, m1_w = GaussQuadrature.legendre(FT,qm1,GaussQuadrature.both)
   m1_phi, m1_phir = jacobip(alpha, beta, deg, m1_r)
   trans = inv(m1_phi)
   return new{T,FT}(qm1, m1_r, m1_w, m1_phi, trans)
end # function B1d / inner constructor
#-----------------------------------------
end # struct B1d
#---------------------------------------------
"""
    jacobip(alpha::T, beta::T, np::T, x::Vector{FT}; dflag=true) where {T <: Integer, FT <: AbstractFloat}

Computes the first np+1 jacobi polynomials and their derivatives on a 1D grid x. The returned matrix sizes are (nx x np+1).
Computation of derivatives can be suppressed using the dflag.  
"""
function jacobip(alpha::T, beta::T, np::T, x::Vector{FT}; dflag=true) where {T <: Integer, FT <: AbstractFloat}
    nx = length(x)
    a  = Vector{FT}(undef,4)
    V  = Array{FT}(undef, nx, np+1)

    if dflag
        DV = Array{FT}(undef, nx, np+1)
    end
    #-------------------------------------------------------
    @assert np ≥ 0
    V  .= 0.0
    if dflag 
        DV .= 0.0
    end
    
    V[:,1] .= 1.0

    if np > 0
        V[:,2] .= 0.5 .* ( alpha .- beta .+ (alpha .+ beta .+ 2.0) .* x )
        if dflag
            DV[:,2] .= 0.5 .* ( alpha .+ beta .+ 2.0 )
        end

	    if np>1
            for i in 2:np 
                a[1] = (2*i)*(i+alpha+beta)*(2*i + alpha + beta -2)
    		    a[2] = (2*i + alpha + beta - 1)*(alpha*alpha - beta*beta)						 
	    	    a[3] = (2*i + alpha + beta - 2)*(2*i + alpha + beta - 1)*(2*i + alpha + beta)  
		        a[4] = 2*(i+alpha-1)*(i+beta-1)*(2*i + alpha + beta)	

                V[:,i+1] .= ( (a[2] .+ a[3].*x) .*  V[:,i] .- a[4].*V[:,i-1] ) ./ a[1]
                if dflag
                    DV[:,i+1] .= ( (a[2] .+ a[3].*x) .* DV[:,i] .+ a[3].*V[:,i] - a[4] .* DV[:,i-1] ) ./ a[1]
                end # for if dflag
            end # for i
        end # if np > 1
    end # if np > 0
    #---------------------------------------------------------
    if dflag 
        return V, DV
    else 
        return V
    end
end # jacobip
#---------------------------------------------
end # module B1d
