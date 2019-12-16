module Interpolation

using MPI
import GaussQuadrature
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Mesh.Geometry
using CLIMA.Mesh.Elements
using CLIMA.Mesh.Tens
using StaticArrays

export Interpolation_Brick, interpolate_brick!

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
end # module interploation
