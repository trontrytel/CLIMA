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

    xmin::FT;  ymin::FT;  zmin::FT # domain bounds, min
    xmax::FT;  ymax::FT;  zmax::FT # & max
    xres::FT;  yres::FT;  zres::FT # respective resolutions for the uniform grid

    x::Vector{Vector{FT}} # unique x coordinates of interpolation points within each element
    y::Vector{Vector{FT}} # unique y coordinates of interpolation points within each element
    z::Vector{Vector{FT}} # unique z coordinates of interpolation points within each element

    ξ1::Vector{Vector{FT}} # unique ξ1 coordinates of interpolation points within each element 
    ξ2::Vector{Vector{FT}} # unique ξ2 coordinates of interpolation points within each element 
    ξ3::Vector{Vector{FT}} # unique ξ3 coordinates of interpolation points within each element 
  
    xg::Vector{Vector{FT}} # x coordinates of portion of interpolation grid embedded within each element
    yg::Vector{Vector{FT}} # y coordinates of portion of interpolation grid embedded within each element
    zg::Vector{Vector{FT}} # z coordinates of portion of interpolation grid embedded within each element

    V::Vector{Vector{FT}}  # interpolated variable within each element

function Interpolation_Brick(grid::DiscontinuousSpectralElementGrid, xres::FT, yres::FT, zres::FT) where FT <: AbstractFloat

    xmin, xmax = FT(minimum( grid.topology.elemtocoord[1,:,:] )), FT(maximum( grid.topology.elemtocoord[1,:,:] ))
    ymin, ymax = FT(minimum( grid.topology.elemtocoord[2,:,:] )), FT(maximum( grid.topology.elemtocoord[2,:,:] ))
    zmin, zmax = FT(minimum( grid.topology.elemtocoord[3,:,:] )), FT(maximum( grid.topology.elemtocoord[3,:,:] ))

    xgrd, ygrd, zgrd = range(xmin, xmax, step=xres), range(ymin, ymax, step=yres), range(zmin, zmax, step=zres) 
    #-----------------------------------------------------------------------------------
    El = grid.topology.realelems # Element (numbers) on the local processor

    np, xst, xen, yst, yen, zst, zen = ntuple( i -> zeros(Integer, length(El)), 7)
    nx, ny, nz = ntuple( i -> zeros(Integer, length(El)), 3)
    xminl, xmaxl, yminl, ymaxl, zminl, zmaxl = ntuple( i -> zeros(FT, length(El)), 6) # x,y,z limits for each brick element
    #-----------------------------------------------------------------------------------
    ectr = 1
    for el in El
      xminl[ectr], xmaxl[ectr] = minimum( grid.topology.elemtocoord[1,:,ectr] ), maximum( grid.topology.elemtocoord[1,:,ectr] )
      yminl[ectr], ymaxl[ectr] = minimum( grid.topology.elemtocoord[2,:,ectr] ), maximum( grid.topology.elemtocoord[2,:,ectr] )
      zminl[ectr], zmaxl[ectr] = minimum( grid.topology.elemtocoord[3,:,ectr] ), maximum( grid.topology.elemtocoord[3,:,ectr] )

      xst[ectr], xen[ectr] = findfirst( temp -> temp ≥ xminl[ectr], xgrd ), findlast( temp -> temp ≤ xmaxl[ectr], xgrd )
      yst[ectr], yen[ectr] = findfirst( temp -> temp ≥ yminl[ectr], ygrd ), findlast( temp -> temp ≤ ymaxl[ectr], ygrd )
      zst[ectr], zen[ectr] = findfirst( temp -> temp ≥ zminl[ectr], zgrd ), findlast( temp -> temp ≤ zmaxl[ectr], zgrd )

      nx[ectr], ny[ectr], nz[ectr] = xen[ectr]-xst[ectr]+1, yen[ectr]-yst[ectr]+1, zen[ectr]-zst[ectr]+1

      np[ectr] = nx[ectr]*ny[ectr]*nz[ectr]

      ectr += 1 
    end # el loop
    #-----------------------------------------------------------------------------------
    x = map( i -> zeros(FT,i), nx)
    y = map( i -> zeros(FT,i), ny)
    z = map( i -> zeros(FT,i), nz)
 
    ξ1 = map( i -> zeros(FT,i), nx)
    ξ2 = map( i -> zeros(FT,i), ny)
    ξ3 = map( i -> zeros(FT,i), nz)

    xg = map( i -> zeros(FT,i), np)
    yg = map( i -> zeros(FT,i), np)
    zg = map( i -> zeros(FT,i), np)

    V = map( i -> zeros(FT,i), np)
    #-----------------------------------------------------------------------------------
    ectr = 1
    for el in El
        [ x[ectr][i] = xgrd[ xst[ectr] + i - 1] for i in 1:nx[ectr] ]
        [ y[ectr][i] = ygrd[ yst[ectr] + i - 1] for i in 1:ny[ectr] ]
        [ z[ectr][i] = zgrd[ zst[ectr] + i - 1] for i in 1:nz[ectr] ]

        [ ξ1[ectr][i] = 2.0 * ( x[ectr][i] - xminl[ectr] ) / (xmaxl[ectr]-xminl[ectr]) -  1.0 for i in 1:nx[ectr] ]
        [ ξ2[ectr][i] = 2.0 * ( y[ectr][i] - yminl[ectr] ) / (ymaxl[ectr]-yminl[ectr]) -  1.0 for i in 1:ny[ectr] ]
        [ ξ3[ectr][i] = 2.0 * ( z[ectr][i] - zminl[ectr] ) / (zmaxl[ectr]-zminl[ectr]) -  1.0 for i in 1:nz[ectr] ]

        ctr = 1

        for k in 1:nz[ectr]
            for j in 1:ny[ectr]
                for i in 1:nx[ectr]
                    xg[ectr][ctr] = x[ectr][i]
                    yg[ectr][ctr] = y[ectr][j]
                    zg[ectr][ctr] = z[ectr][k]
                    ctr += 1
                end
            end
        end
      ectr += 1
    end # el loop
    #-----------------------------------------------------------------------------------
    Nel = length(El)
    return new{FT}(El, Nel, xmin, ymin, zmin, xmax, ymax, zmax, xres, yres, zres, x, y, z, ξ1, ξ2, ξ3, xg, yg, zg, V)

end # function Interpolation_Brick -> inner constructor
#--------------------------------------------------------
end # struct Interpolation_Brick 
#--------------------------------------------------------
function interpolate_brick!(intrp_brck::Interpolation_Brick, sv::AbstractArray{FT}, st_no::T, poly_order::T) where {T <: Integer, FT <: AbstractFloat}

  qm1 = poly_order + 1
  m1_r, m1_w = GaussQuadrature.legendre(FT,qm1,GaussQuadrature.both)

  #-----for each element elno 
  for el in 1:intrp_brck.Nel

    lag = @view sv[:,st_no,el]

    g_phir = Elements.interpolationmatrix(m1_r, intrp_brck.ξ1[el])
    g_phis = Elements.interpolationmatrix(m1_r, intrp_brck.ξ2[el])
    g_phit = Elements.interpolationmatrix(m1_r, intrp_brck.ξ3[el])

    tenspxv_hex!(g_phir, g_phis, g_phit, false, lag, intrp_brck.V[el]) 
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

    lat_grd::Array{FT}                        # lat grid locations
    long_grd::Array{FT}                       # long grid locations
    rad_grd::Array{FT}                        # rad grid locations

    n_lat::T; n_long::T; n_rad::T;            # # of lat, long & rad grid locations

    el_grd                          # element containing the grid point


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
    if ( ctr > max_it )
      error("invert_trilinear_mapping_hex: Newton-Raphson not converging to desired tolerance after max_it iterations")
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
