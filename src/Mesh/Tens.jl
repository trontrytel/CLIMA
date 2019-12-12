module Tens
using LinearAlgebra
export tenspxv_hex!

"""
    tenspxv_hex!(phir::Array{FT,2}, phis::Array{FT,2}, phit::Array{FT,2},
                     tr_opt::Bool, vin::AbstractArray{FT}, vout::AbstractArray{FT} ) where FT <: AbstractFloat 

Computes the tensor product x vector voutʳˢᵗ = phirʳⁱ phisˢʲ phitᵗᵏ vinⁱʲᵏ	 or its transpose version. 
"""
function tenspxv_hex!(phir::Array{FT,2}, phis::Array{FT,2}, phit::Array{FT,2},
                     tr_opt::Bool,                                               #transpose (true/false) 
                     vin::AbstractArray{FT}, vout::AbstractArray{FT} ) where FT <: AbstractFloat #input and output vectors 
                     

  si, sj, sk = size(phir,2), size(phis,2), size(phit,2)
  sr, ss, st = size(phir,1), size(phis,1), size(phit,1)
  alpha = FT(1.0)
  beta  = FT(0.0)  

  if tr_opt # transpose
    #----------------------------------------------------------------
    h1 = zeros(FT, ss*st, si)  # this is a test version
    h2 = zeros(FT, st*si, sj)  # memory allocations within the function will be removed
    #----------------------------------------------------------------
    LinearAlgebra.BLAS.gemm!('T', 'N', alpha, reshape(vin,sr, ss*st), phir, beta, h1)
    LinearAlgebra.BLAS.gemm!('T', 'N', alpha, reshape(h1, ss, st*si), phis, beta, h2)
    LinearAlgebra.BLAS.gemm!('T', 'N', alpha, reshape(h2, st, si*sj), phit, beta, reshape(vout, si*sj, sk))
    #---------------------------------------------------------------
  else # no transpose
    #---------------------------------------------------------------
    h1 = zeros(FT, sj*sk, sr)  # this is for testing only
    h2 = zeros(FT, sk*sr, ss)  # memory allocations within the function will be removed
    #---------------------------------------------------------------
    LinearAlgebra.BLAS.gemm!('T', 'T', alpha, reshape(vin,si, sj*sk), phir, beta, h1)
    LinearAlgebra.BLAS.gemm!('T', 'T', alpha, reshape(h1, sj, sk*sr), phis, beta, h2)
    LinearAlgebra.BLAS.gemm!('T', 'T', alpha, reshape(h2, sk, sr*ss), phit, beta, reshape(vout, sr*ss, st))
    #---------------------------------------------------------------
  end

end # function tenspxv_hex!
#-----------------------------
end # module Tens
