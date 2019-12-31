module Interpolation

using MPI
import GaussQuadrature
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Mesh.Geometry
using CLIMA.Mesh.Elements
using CLIMA.Mesh.Tens
using LinearAlgebra
using StaticArrays

export Interpolation_Brick, interpolate_brick!, Interpolation_Cubed_Sphere, invert_trilear_mapping_hex

#--------------------------------------------------------

# This interpolation routine works for a brick, where stretching/compression happens only along the x, y & z axis.
# Here x = X(ξ1); y = Y(ξ2); z = Z(ξ3)

struct Interpolation_Brick{FT <:AbstractFloat}
    El::UnitRange{Int64}
    Nel::Integer

    xbnd::Array{FT}        # domain bounds, [2(min/max),ndim]
    xres::Array{FT}        # resolutions in x1, x2, x3 directions for the uniform grid

    ξ1::Vector{Vector{FT}} # unique ξ1 coordinates of interpolation points within each element 
    ξ2::Vector{Vector{FT}} # unique ξ2 coordinates of interpolation points within each element 
    ξ3::Vector{Vector{FT}} # unique ξ3 coordinates of interpolation points within each element 
  
    x::Vector{Array{FT}}   # x[elno][3,npts] -> (x1,x2,x3) coordinates of portion of interpolation grid embedded within each element 

    V::Vector{Vector{FT}}  # interpolated variable within each element
#--------------------------------------------------------
function Interpolation_Brick(grid::DiscontinuousSpectralElementGrid, xres::Array{FT}) where FT <: AbstractFloat
    T = Integer
    ndim = 3
    xbnd = zeros(FT, 2, ndim) # domain bounds (min,max) in each dimension
    for dim in 1:ndim 
        xbnd[1,dim], xbnd[2,dim] = FT(minimum( grid.topology.elemtocoord[dim,:,:] )), FT(maximum( grid.topology.elemtocoord[dim,:,:] )) 
    end
    x1g, x2g, x3g = range(xbnd[1,1], xbnd[2,1], step=xres[1]), range(xbnd[1,2], xbnd[2,2], step=xres[2]), range(xbnd[1,3], xbnd[2,3], step=xres[3]) 
    #-----------------------------------------------------------------------------------
    El =  grid.topology.realelems # Element (numbers) on the local processor
    Nel = length(El)
    n123  = zeros(T,    ndim)     # # of unique ξ1, ξ2, ξ3 points in each cell
    xsten = zeros(T, 2, ndim)     # x1, x2, x3 start and end for each brick element
    xbndl = zeros(T, 2, ndim)     # location of x1,x2,x3 limits (min,max) for each brick element

    ξ1, ξ2, ξ3 = map( i -> zeros(FT,i), zeros(T,Nel)), map( i -> zeros(FT,i), zeros(T,Nel)), map( i -> zeros(FT,i), zeros(T,Nel))
    V  = map( i -> zeros(FT,i), zeros(T,Nel))
    x  = map( i -> zeros(FT,ndim,i), zeros(T,Nel)) # interpolation grid points embedded in each cell 
    #-----------------------------------------------------------------------------------
    for el in 1:Nel
        for (ξ,xg,dim) in zip((ξ1,ξ2,ξ3), (x1g, x2g, x3g), 1:ndim)
            xbndl[1,dim], xbndl[2,dim] = minimum( grid.topology.elemtocoord[dim,:,el] ), maximum( grid.topology.elemtocoord[dim,:,el] )
            xsten[1,dim], xsten[2,dim] = findfirst( temp -> temp ≥ xbndl[1,dim], xg ), findlast( temp -> temp ≤ xbndl[2,dim], xg )
            n123[dim] = xsten[2,dim] - xsten[1,dim]
            ξ[el] = Array{FT}(undef,n123[dim])
            [ ξ[el][i] = 2.0 * ( xg[ xsten[1,dim] + i - 1] - xbndl[1,dim] ) / (xbndl[2,dim]-xbndl[1,dim]) -  1.0 for i in 1:n123[dim] ]
        end

        x[el] = zeros(FT,ndim,prod(n123))
        V[el] = zeros(FT,prod(n123))

        ctr = 1

        for k in 1:n123[3], j in 1:n123[2], i in 1:n123[1]
            x[el][1,ctr] = x1g[ xsten[1,1] + i - 1 ]
            x[el][2,ctr] = x2g[ xsten[1,2] + j - 1 ]
            x[el][3,ctr] = x3g[ xsten[1,3] + k - 1 ]
            ctr += 1
        end
    end # el loop
    #-----------------------------------------------------------------------------------
    return new{FT}(El, Nel, xbnd, xres, ξ1, ξ2, ξ3, x, V)

end # function Interpolation_Brick -> inner constructor
#--------------------------------------------------------
end # struct Interpolation_Brick 
#--------------------------------------------------------
function interpolate_brick!(intrp_brck::Interpolation_Brick, sv::AbstractArray{FT}, st_no::T, poly_order::T) where {T <: Integer, FT <: AbstractFloat}

    qm1 = poly_order + 1
    m1_r, m1_w = GaussQuadrature.legendre(FT,qm1,GaussQuadrature.both)
    #-----for each element elno 
    for el in 1:intrp_brck.Nel
        if length(intrp_brck.ξ1[el]) > 0
            lag    = @view sv[:,st_no,el]
            g_phir = Elements.interpolationmatrix(m1_r, intrp_brck.ξ1[el])
            g_phis = Elements.interpolationmatrix(m1_r, intrp_brck.ξ2[el])
            g_phit = Elements.interpolationmatrix(m1_r, intrp_brck.ξ3[el])

            tenspxv_hex!(g_phir, g_phis, g_phit, false, lag, intrp_brck.V[el]) 
        end
  #--------------------
  end
end
#--------------------------------------------------------
#--------------------------------------------------------
struct Interpolation_Cubed_Sphere{T <: Integer, FT <: AbstractFloat}

    El::UnitRange{Int64}

    Nel::T

    lat_min::FT;  long_min::FT;  rad_min::FT; # domain bounds, min
    lat_max::FT;  long_max::FT;  rad_max::FT; # domain bounds, max
    lat_res::FT;  long_res::FT;  rad_res::FT; # respective resolutions for the uniform grid

#    x::Vector{Vector{FT}} # unique x coordinates of interpolation points within each element
#    y::Vector{Vector{FT}} # unique y coordinates of interpolation points within each element
#    z::Vector{Vector{FT}} # unique z coordinates of interpolation points within each element

#    ξ1::Vector{Vector{FT}} # unique ξ1 coordinates of interpolation points within each element 
#    ξ2::Vector{Vector{FT}} # unique ξ2 coordinates of interpolation points within each element 
#    ξ3::Vector{Vector{FT}} # unique ξ3 coordinates of interpolation points within each element 


    lat_grd::Array{FT}                        # lat grid locations
    long_grd::Array{FT}                       # long grid locations
    rad_grd::Array{FT}                        # rad grid locations

    n_lat::T; n_long::T; n_rad::T;            # # of lat, long & rad grid locations

    el_grd                                    # element containing the grid point
  #--------------------------------------------------------
  function Interpolation_Cubed_Sphere(grid::DiscontinuousSpectralElementGrid, vert_range::AbstractArray{FT}, lat_res::FT, long_res::FT, rad_res::FT) where {FT <: AbstractFloat}
    T = Integer

    toler1 = eps(FT) * vert_range[1] * 2.0 # tolerance for unwarp function
    toler2 = eps(FT) * 4.0                 # tolerance 

    lat_min,   lat_max = FT(0.0), FT(π)                 # inclination/zeinth angle range
    long_min, long_max = FT(0.0), FT(2*π)  			    # azimuthal angle range
    rad_min,   rad_max = vert_range[1], vert_range[end] # radius range

    El = grid.topology.realelems # Element (numbers) on the local processor
    Nel = length(El)

    Nel_global = length(grid.topology.elems)
    nvert = length(vert_range) - 1              # # of elements in vertical direction
    nhor  = T( sqrt( Nel_global / nvert / 6) )  # # of elements in horizontal direction
    nblck = nhor * nhor * nvert
    Δh = 2.0 / nhor                             # horizontal grid spacing in unwarped grid

    lat_grd, long_grd, rad_grd = range(lat_min, lat_max, step=lat_res), range(long_min, long_max, step=long_res), range(rad_min, rad_max, step=rad_res) 

    n_lat, n_long, n_rad = T(length(lat_grd)), T(length(long_grd)), T(length(rad_grd))

    el_grd = Array{T}(undef, n_rad, n_lat, n_long) # element containing the lat/long/rad grid point 
    uw_grd = zeros(FT, 3, 1)
    #---------------------------------------------- 

    flip_ord = invperm(grid.topology.origsendorder) # to account for reordering of elements after the partitioning process 

    reorder = zeros(T, 6*nvert*nhor*nhor)

    for i in 1:length(flip_ord), j in 1:nvert
        reorder[j + (i-1)*nvert] = (flip_ord[i] - 1)*nvert + j
    end

    temp = invperm(grid.topology.origsendorder) 

    for k in 1:n_long
      for j in 1:n_lat
        for i in 1:n_rad 
          x1_grd = rad_grd[i] * sin(lat_grd[j]) * cos(long_grd[k]) # inclination -> latitude; azimuthal -> longitude.
          x2_grd = rad_grd[i] * sin(lat_grd[j]) * sin(long_grd[k]) # inclination -> latitude; azimuthal -> longitude.
          x3_grd = rad_grd[i] * cos(lat_grd[j])   

          uw_grd[1], uw_grd[2], uw_grd[3] = Topologies.cubedshellunwarp(x1_grd, x2_grd, x3_grd) # unwarping from sphere to cubed shell
          rad = FT(maximum(abs.(uw_grd)))
          #--------------------------------
          x1_uw2_grd = uw_grd[1] / rad # unwrapping cubed shell on to a 2D grid (in 3D space, -1 to 1 cube)
          x2_uw2_grd = uw_grd[2] / rad
          x3_uw2_grd = uw_grd[3] / rad
          #--------------------------------
          if rad ≤ vert_range[1]       # accounting for minor rounding errors from unwarp function at boundaries 
            vert_range[1] - rad < toler1 ? l_nrm = 1 :  error("fatal error, rad lower than inner radius")
          elseif rad ≥ vert_range[end] # accounting for minor rounding errors from unwarp function at boundaries 
            rad - vert_range[end] < toler1 ? l_nrm = nvert : error("fatal error, rad greater than outer radius")
          else                         # normal scenario
            l_nrm = findfirst( X -> X .- rad .> 0.0, vert_range ) - 1 # identify stack bin 
          end
          #--------------------------------
          if     abs(x1_uw2_grd + 1) < toler2 # face 1 (x1 == -1 plane)
		    l2 = min(div(x2_uw2_grd + 1, Δh) + 1, nhor)
		    l3 = min(div(x3_uw2_grd + 1, Δh) + 1, nhor)
            el_grd[i,j,k] = reorder[T(l_nrm + (nhor-l2)*nvert + (l3-1)*nvert*nhor)] 
          elseif abs(x2_uw2_grd + 1) < toler2 # face 2 (x2 == -1 plane)
		    l1 = min(div(x1_uw2_grd + 1, Δh) + 1, nhor)
		    l3 = min(div(x3_uw2_grd + 1, Δh) + 1, nhor)
            el_grd[i,j,k] = reorder[T(l_nrm + (l1-1)*nvert + (l3-1)*nvert*nhor + nblck*1)]
          elseif abs(x1_uw2_grd - 1) < toler2 # face 3 (x1 == +1 plane)
		    l2 = min(div(x2_uw2_grd + 1, Δh) + 1, nhor)
		    l3 = min(div(x3_uw2_grd + 1, Δh) + 1, nhor)
            el_grd[i,j,k] = reorder[T(l_nrm + (l2-1)*nvert + (l3-1)*nvert*nhor + nblck*2 )]
          elseif abs(x3_uw2_grd - 1) < toler2 # face 4 (x3 == +1 plane)
  		    l1 = min(div(x1_uw2_grd + 1, Δh) + 1, nhor)
		    l2 = min(div(x2_uw2_grd + 1, Δh) + 1, nhor)
            el_grd[i,j,k] = reorder[T(l_nrm + (l1-1)*nvert + (l2-1)*nvert*nhor + nblck*3)]
          elseif abs(x2_uw2_grd - 1) < toler2 # face 5 (x2 == +1 plane)
	        l1 = min(div(x1_uw2_grd + 1, Δh) + 1, nhor)
		    l3 = min(div(x3_uw2_grd + 1, Δh) + 1, nhor)
            el_grd[i,j,k] = reorder[T(l_nrm + (l1-1)*nvert + (nhor-l3)*nvert*nhor + nblck*4 )]
          elseif abs(x3_uw2_grd + 1) < toler2 # face 6 (x3 == -1 plane)
		    l1 = min(div(x1_uw2_grd + 1, Δh) + 1, nhor)
		    l2 = min(div(x2_uw2_grd + 1, Δh) + 1, nhor)
            el_grd[i,j,k] = reorder[T(l_nrm + (l1-1)*nvert + (nhor-l2)*nvert*nhor + nblck*5)]
          else
            error("error: unwrapped grid does on lie on any of the 6 faces")
          end
          #--------------------------------
          #---testing-----------------------------
#          it = 1; jt = 1; kt = 1
#          if i==it && j==jt && k == kt
#            println("testing")
            elm = el_grd[i,j,k]
            X = grid.topology.elemtocoord[1,:,el_grd[i,j,k]]
            Y = grid.topology.elemtocoord[2,:,el_grd[i,j,k]]
            Z = grid.topology.elemtocoord[3,:,el_grd[i,j,k]]
            ξ = invert_trilear_mapping_hex(X, Y, Z, uw_grd)
            println("Elem # ",el_grd[i,j,k],"; ξ = ", ξ)
#          end
          #--------------------------------
      end
    end
  end
 

  return new{T, FT}(El, Nel, lat_min, long_min, rad_min, lat_max, long_max, rad_max, lat_res, long_res, rad_res, 
                    lat_grd, long_grd, rad_grd, n_lat, n_long, n_rad, el_grd)
    #-----------------------------------------------------------------------------------
end # Inner constructor function Interpolation_Cubed_Sphere
end # structure Interpolation_Cubed_Sphere
#--------------------------------------------------------
# This function computes (ξ1,ξ2,ξ3) given (x1,x2,x3) and the (8) vertex coordinates of a Hexahedron
function invert_trilear_mapping_hex(X1::Array{FT}, X2::Array{FT}, X3::Array{FT}, x::Array{FT}) where FT <: AbstractFloat 
  tol    = eps(FT) * 4.0 # tolerance for Newton-Raphson solver
  max_it = 10            # maximum # of iterations
  ξ      = zeros(FT,3,1) # initial guess => cell centroid

  d   = trilinear_map(ξ, X1, X2, X3) - x
  err = norm(d)
  ctr = 0 
  #---Newton-Raphson iterations---------------------------
  while err > tol
    trilinear_map_IJac_x_vec!(ξ, X1, X2, X3, d)
    ξ .-= d
    d = trilinear_map(ξ, X1, X2, X3) - x
    err = norm(d)
    ctr += 1
    if ctr > max_it
      error("invert_trilinear_mapping_hex: Newton-Raphson not converging to desired tolerance after max_it = ", max_it," iterations")
    end
  end
  #-------------------------------------------------------
  return ξ
end
#--------------------------------------------------------
function trilinear_map(ξ::Array{FT}, x1v::Array{FT}, x2v::Array{FT}, x3v::Array{FT}) where FT <: AbstractFloat
  x = Array{FT}(undef,3)
  for (vert,dim) = zip((x1v,x2v,x3v),1:3)
    x[dim] = ((1 - ξ[1]) * (1 - ξ[2]) * (1 - ξ[3]) * vert[1] + (1 + ξ[1]) * (1 - ξ[2]) * (1 - ξ[3]) * vert[2] +
              (1 - ξ[1]) * (1 + ξ[2]) * (1 - ξ[3]) * vert[3] + (1 + ξ[1]) * (1 + ξ[2]) * (1 - ξ[3]) * vert[4] +
              (1 - ξ[1]) * (1 - ξ[2]) * (1 + ξ[3]) * vert[5] + (1 + ξ[1]) * (1 - ξ[2]) * (1 + ξ[3]) * vert[6] +
              (1 - ξ[1]) * (1 + ξ[2]) * (1 + ξ[3]) * vert[7] + (1 + ξ[1]) * (1 + ξ[2]) * (1 + ξ[3]) * vert[8] )/ 8.0
  end
  return x
end
#--------------------------------------------------------
function trilinear_map_IJac_x_vec!(ξ::Array{FT}, x1v::Array{FT}, x2v::Array{FT}, x3v::Array{FT}, v::Array{FT}) where FT <: AbstractFloat
  Jac = zeros(FT, 3,3)
  for (vert,dim) = zip((x1v,x2v,x3v),1:3)
    Jac[dim,1] = ((-1) * (1 - ξ[2]) * (1 - ξ[3]) * vert[1] + ( 1) * (1 - ξ[2]) * (1 - ξ[3]) * vert[2] +
                  (-1) * (1 + ξ[2]) * (1 - ξ[3]) * vert[3] + (+1) * (1 + ξ[2]) * (1 - ξ[3]) * vert[4] +
                  (-1) * (1 - ξ[2]) * (1 + ξ[3]) * vert[5] + (+1) * (1 - ξ[2]) * (1 + ξ[3]) * vert[6] +
                  (-1) * (1 + ξ[2]) * (1 + ξ[3]) * vert[7] + (+1) * (1 + ξ[2]) * (1 + ξ[3]) * vert[8] )/ 8.0

    Jac[dim,2] = ((1 - ξ[1]) * (-1) * (1 - ξ[3]) * vert[1] + (1 + ξ[1]) * (-1) * (1 - ξ[3]) * vert[2] +
                  (1 - ξ[1]) * (+1) * (1 - ξ[3]) * vert[3] + (1 + ξ[1]) * (+1) * (1 - ξ[3]) * vert[4] +
                  (1 - ξ[1]) * (-1) * (1 + ξ[3]) * vert[5] + (1 + ξ[1]) * (-1) * (1 + ξ[3]) * vert[6] +
                  (1 - ξ[1]) * (+1) * (1 + ξ[3]) * vert[7] + (1 + ξ[1]) * (+1) * (1 + ξ[3]) * vert[8] )/ 8.0

    Jac[dim,3] = ((1 - ξ[1]) * (1 - ξ[2]) * (-1) * vert[1] + (1 + ξ[1]) * (1 - ξ[2]) * (-1) * vert[2] +
                  (1 - ξ[1]) * (1 + ξ[2]) * (-1) * vert[3] + (1 + ξ[1]) * (1 + ξ[2]) * (-1) * vert[4] +
                  (1 - ξ[1]) * (1 - ξ[2]) * (+1) * vert[5] + (1 + ξ[1]) * (1 - ξ[2]) * (+1) * vert[6] +
                  (1 - ξ[1]) * (1 + ξ[2]) * (+1) * vert[7] + (1 + ξ[1]) * (1 + ξ[2]) * (+1) * vert[8] )/ 8.0
  end
  LinearAlgebra.LAPACK.gesv!(Jac,v)
  return nothing 
end
#--------------------------------------------------------
end # module interploation
