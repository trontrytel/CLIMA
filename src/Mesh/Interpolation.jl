module Interpolation
using CLIMA
using MPI
import GaussQuadrature
using CLIMA.Mesh.Topologies
using CLIMA.Mesh.Grids
using CLIMA.Mesh.Geometry
using CLIMA.Mesh.Elements
using LinearAlgebra
using StaticArrays

#-------------------
#using Requires
#@init @require CUDAnative = "be33ccc6-a3ff-5ff2-a52e-74243cff1e17" begin
# using .CUDAnative
#end

using CUDAdrv
using CUDAnative
using CuArrays
using CUDAapi
#-------------------
using GPUifyLoops

export InterpolationBrick, interpolate_brick!, InterpolationCubedSphere, invert_trilear_mapping_hex, interpolate_cubed_sphere!, interpolate_cubed_sphere_CUDA!


#interpolate_brick_v2!
#--------------------------------------------------------
"""
    InterpolationBrick(grid::DiscontinuousSpectralElementGrid, xres, ::Type{FT}) where FT <: AbstractFloat

This interpolation structure and the corresponding functions works for a brick, where stretching/compression happens only along the x, y & z axis.
Here x1 = X1(ξ1); x2 = X2(ξ2); x3 = X3(ξ3)

# input for the inner constructor
 - `grid` DiscontinousSpectralElementGrid
 - `xres` Resolution of the interpolation grid in x1, x2 and x3 directions
"""
struct InterpolationBrick{FT <:AbstractFloat}
    realelems::UnitRange{Int64}
    poly_order::Int

    xbnd::Array{FT,2}      # domain bounds, [2(min/max),ndim]
    xres::Array{FT,1}      # resolutions in x1, x2, x3 directions for the uniform grid

    #ξ1::Vector{Vector{FT}} # unique ξ1 coordinates of interpolation points within each element 
    #ξ2::Vector{Vector{FT}} # unique ξ2 coordinates of interpolation points within each element 
    #ξ3::Vector{Vector{FT}} # unique ξ3 coordinates of interpolation points within each element 

    ξ1::AbstractVector{FT} # unique ξ1 coordinates of interpolation points within each element 
    ξ2::AbstractVector{FT} # unique ξ2 coordinates of interpolation points within each element 
    ξ3::AbstractVector{FT} # unique ξ3 coordinates of interpolation points within each element 

    o1::AbstractVector{Int} # offsets for each element for ξ1
    o2::AbstractVector{Int}
    o3::AbstractVector{Int}

    flg1::AbstractVector{Int} # flags for ξ1
    flg2::AbstractVector{Int}
    flg3::AbstractVector{Int}

    fac1::AbstractVector{FT} # factor for ξ1 
    fac2::AbstractVector{FT}
    fac3::AbstractVector{FT}
  
    x::Vector{Array{FT,2}} # x[elno][3,npts] -> (x1,x2,x3) coordinates of portion of interpolation grid embedded within each element 

    v::AbstractArray{FT,1} # interpolated variable within each element
    offset::AbstractVector{Int}    # offsets for each element for v

    m1_r::AbstractVector{FT}   # GLL points
    m1_w::AbstractVector{FT}   # GLL weights
    wb::AbstractVector{FT}     # Barycentric weights

#--------------------------------------------------------
    function InterpolationBrick(grid::DiscontinuousSpectralElementGrid{FT}, xres) where FT <: AbstractFloat
		DA = CLIMA.array_type()                    # device array
        T = Int
        poly_order = polynomialorder(grid)
        qm1 = poly_order + 1
        ndim = 3
        toler = FT(eps(FT) * 4.0)                 # tolerance 
        xbnd = zeros(FT, 2, ndim) # domain bounds (min,max) in each dimension
        for dim in 1:ndim 
            xbnd[1,dim], xbnd[2,dim] = extrema(grid.topology.elemtocoord[dim,:,:])
        end
        x1g = range(xbnd[1,1], xbnd[2,1], step=xres[1])
        x2g = range(xbnd[1,2], xbnd[2,2], step=xres[2])
        x3g = range(xbnd[1,3], xbnd[2,3], step=xres[3]) 
        #-----------------------------------------------------------------------------------
        realelems = grid.topology.realelems  # Element (numbers) on the local processor
        Nel       = length(realelems)
        offset    = Vector{Int}(undef,Nel+1) # offsets for the interpolated variable
        o1_d      = Vector{Int}(undef,Nel+1) # offsets for ξ1
        o2_d      = Vector{Int}(undef,Nel+1) # offsets for ξ2
        o3_d      = Vector{Int}(undef,Nel+1) # offsets for ξ3
        n123      = zeros(T,    ndim)        # # of unique ξ1, ξ2, ξ3 points in each cell
        xsten     = zeros(T, 2, ndim)        # x1, x2, x3 start and end for each brick element
        xbndl     = zeros(FT,2, ndim)        # x1,x2,x3 limits (min,max) for each brick element

        ξ1 = map( i -> zeros(FT,i), zeros(T,Nel))
        ξ2 = map( i -> zeros(FT,i), zeros(T,Nel))
        ξ3 = map( i -> zeros(FT,i), zeros(T,Nel))

        x  = map( i -> zeros(FT,ndim,i), zeros(T,Nel)) # interpolation grid points embedded in each cell 

        offset[1] = 0; o1_d[1] = 0; o2_d[1] = 0; o3_d[1] = 0;

        #-----------------------------------------------------------------------------------
        for el in 1:Nel
            for (ξ,xg,dim,offv) in zip((ξ1,ξ2,ξ3), (x1g, x2g, x3g), 1:ndim,(o1_d,o2_d,o3_d))
                xbndl[1,dim], xbndl[2,dim] = extrema(grid.topology.elemtocoord[dim,:,el])

                
                st = findfirst( map( (av,bv) -> av && bv, (xg .≥ xbndl[1,dim]), (xg .≤ xbndl[2,dim])) )
           
                if st ≠ nothing
                    xsten[1,dim] = st
                    xsten[2,dim] = findlast( temp -> temp ≤ xbndl[2,dim], xg )
                    n123[dim] = xsten[2,dim] - xsten[1,dim] + 1
                end
                ξ[el] = [ 2 * ( xg[ xsten[1,dim] + i - 1] - xbndl[1,dim] ) / (xbndl[2,dim]-xbndl[1,dim]) -  1 for i in 1:n123[dim] ]
                offv[el+1]= offv[el] + n123[dim]
            end

            x_el         = zeros(FT,ndim,prod(n123))
            offset[el+1] = offset[el] + prod(n123) 

            ctr = 1

            for k in 1:n123[3], j in 1:n123[2], i in 1:n123[1]
                x_el[1,ctr] = x1g[ xsten[1,1] + i - 1 ]
                x_el[2,ctr] = x2g[ xsten[1,2] + j - 1 ]
                x_el[3,ctr] = x3g[ xsten[1,3] + k - 1 ]
                ctr += 1
            end
            x[el] = x_el
        end # el loop
        #-----------------------------------------------------------------------------------
        m1_r, m1_w = GaussQuadrature.legendre(FT,qm1,GaussQuadrature.both)
        wb = Elements.baryweights(m1_r)

        ξ1_d   = Array{FT}(undef,o1_d[end]); ξ2_d   = Array{FT}(undef,o2_d[end]); ξ3_d   = Array{FT}(undef,o3_d[end])
        fac1_d = zeros(FT,o1_d[end]); fac2_d = zeros(FT,o2_d[end]); fac3_d = zeros(FT,o3_d[end])
        flg1_d = zeros(T,o1_d[end]);  flg2_d = zeros(T,o2_d[end]);  flg3_d = zeros(T,o3_d[end])

        for i in 1:Nel
            ctr = 1
            for j in o1_d[i]+1:o1_d[i+1]
                ξ1_d[j] = ξ1[i][ctr]
                l1 = findfirst(X -> abs.(X.-ξ1_d[j]) < toler, m1_r)
                l1 == nothing ? [ fac1_d[j] +=  wb[ib] / (ξ1_d[j]-m1_r[ib]) for ib in 1:qm1 ] : (flg1_d[j] = l1; fac1_d[j] = FT(1);)
                fac1_d[j] = FT(1)/fac1_d[j]
                ctr += 1
            end

            ctr = 1
            for j in o2_d[i]+1:o2_d[i+1]
                ξ2_d[j] = ξ2[i][ctr]
                l2 = findfirst(X -> abs.(X.-ξ2_d[j]) < toler, m1_r)
                l2 == nothing ? [ fac2_d[j] +=  wb[ib] / (ξ2_d[j]-m1_r[ib]) for ib in 1:qm1 ] : (flg2_d[j] = l2; fac2_d[j] = FT(1);)
                fac2_d[j] = FT(1)/fac2_d[j]
                ctr += 1
            end

            ctr = 1
            for j in o3_d[i]+1:o3_d[i+1]
                ξ3_d[j] = ξ3[i][ctr]
                l3 = findfirst(X -> abs.(X.-ξ3_d[j]) < toler, m1_r)
                l3 == nothing ? [ fac3_d[j] +=  wb[ib] / (ξ3_d[j]-m1_r[ib]) for ib in 1:qm1 ] : (flg3_d[j] = l3; fac3_d[j] = FT(1);)
                fac3_d[j] = FT(1)/fac3_d[j]
                ctr += 1
            end

        end

        ξ1_d   = DA(ξ1_d);   ξ2_d   = DA(ξ2_d);   ξ3_d   = DA(ξ3_d)
        o1_d   = DA(o1_d);   o2_d   = DA(o2_d);   o3_d   = DA(o3_d)
        flg1_d = DA(flg1_d); flg2_d = DA(flg2_d); flg3_d = DA(flg3_d)
        fac1_d = DA(fac1_d); fac2_d = DA(fac2_d); fac3_d = DA(fac3_d)
        m1_r   = DA(m1_r);   m1_w   = DA(m1_w);   wb     = DA(wb); offset = DA(offset)
        v = DA( Array{FT}(undef,offset[end]) )
        return new{FT}(realelems, poly_order, xbnd, xres, ξ1_d, ξ2_d, ξ3_d, o1_d, o2_d, o3_d, flg1_d, flg2_d, flg3_d, fac1_d, fac2_d, fac3_d, x, v, offset, m1_r, m1_w, wb)

    end
#--------------------------------------------------------
end # struct InterpolationBrick 
#--------------------------------------------------------
"""
    interpolate_brick!(intrp_brck::InterpolationBrick, sv::AbstractArray{FT}, st_idx::T, poly_order::T) where {T <: Integer, FT <: AbstractFloat}

This interpolation function works for a brick, where stretching/compression happens only along the x, y & z axis.
Here x1 = X1(ξ1); x2 = X2(ξ2); x3 = X3(ξ3)

# input
 - `intrp_brck` Initialized InterpolationBrick structure
 - `sv` State vector
 - `st_idx` # of state vector variable to be interpolated
 - `poly_order` polynomial order for the simulation
"""
function interpolate_brick!(offset::AbstractArray{T,1},  m1_r::AbstractArray{FT,1},   wb::AbstractArray{FT,1}, 
                                ξ1::AbstractArray{FT,1},   ξ2::AbstractArray{FT,1},   ξ3::AbstractArray{FT,1}, 
                                o1::AbstractArray{T,1},    o2::AbstractArray{T,1},    o3::AbstractArray{T,1}, 
                              flg1::AbstractArray{T,1},  flg2::AbstractArray{T,1},  flg3::AbstractArray{T,1}, 
                              fac1::AbstractArray{FT,1}, fac2::AbstractArray{FT,1}, fac3::AbstractArray{FT,1}, 
                                 v::AbstractArray{FT,1},   sv::AbstractArray{FT}, st_idx::T) where {T <: Int, FT <: AbstractFloat}
    DA  = CLIMA.array_type()                    # device array
    qm1 = length(m1_r)
    Nel = length(offset) - 1
    
    device = typeof(sv) <: Array ? CPU() : CUDA()
    println("device = $device")

    if device==CPU()
    	vout    = FT(0)
	    vout_ii = FT(0)
    	vout_ij = FT(0)

    	for el in 1:Nel #-----for each element elno 
	        s1 = o1[el]+1; e1 = o1[el+1]; l1 = e1-s1+1;
    	    s2 = o2[el]+1; e2 = o2[el+1]; l2 = e2-s2+1;  
        	s3 = o3[el]+1; e3 = o3[el+1]; l3 = e3-s3+1;

	        if l1 > 0
    	        f1  = view(flg1,s1:e1); fac1l = view(fac1,s1:e1); ξ1l = view(ξ1,s1:e1);
        	    f2  = view(flg2,s2:e2); fac2l = view(fac2,s2:e2); ξ2l = view(ξ2,s2:e2); 
            	f3  = view(flg3,s3:e3); fac3l = view(fac3,s3:e3); ξ3l = view(ξ3,s3:e3); 
	            lag = view(sv,:,st_idx,el)

    	        for k in 1:l3, j in 1:l2, i in 1:l1 # interpolating point-by-point
            	    vout = 0            
                    f3[k] == 0 ? (ikloop = 1:qm1) : (ikloop = f3[k]:f3[k])
                   	for ik in ikloop #1:qm1
                       	v_ij = 0
	                    f2[j] == 0 ? (ijloop = 1:qm1) : (ijloop = f2[j]:f2[j])
    	                for ij in ijloop #1:qm1
        	                #----------------------------------------------------------------------
            	            if f1[i]==0
                                v_ii = 0
                   	            for ii in 1:qm1
                       	            @inbounds v_ii += lag[ii + (ij-1)*qm1 + (ik-1)*qm1*qm1] *  wb[ii] / (ξ1l[i]-m1_r[ii]) #phir[i,ii]
                           	    end # ii loop
	                        else
    	                        @inbounds v_ii = lag[f1[i] + (ij-1)*qm1 + (ik-1)*qm1*qm1]
        	                end
            	            if f2[j]==0
                                @inbounds v_ij += v_ii * wb[ij] / (ξ2l[j]-m1_r[ij])#phis[j,ij]
                  	        else
                       	        v_ij = v_ii
                           	end
	                        #----------------------------------------------------------------------
    	                end # ij loop
        	            if f3[k]==0
            	            @inbounds vout += v_ij * wb[ik] / (ξ3l[k]-m1_r[ik]) #phit[ik]
                        else
                   	        vout = v_ij
                       	end
	                end # ik loop
    	            @inbounds v[ offset[el] + i + (j-1)*l1 + (k-1)*l1*l2] = (vout * fac1l[i] * fac2l[j] * fac3l[k]) 
            	end # i,j,k loop
	        end
	    end
    #--------------------
    else
		@cuda threads=(qm1,qm1) blocks=Nel shmem=qm1*(qm1+2)*sizeof(FT) interpolate_brick_CUDA!(offset,  m1_r, wb, ξ1, ξ2, ξ3, 
                                                      o1, o2, o3, flg1, flg2, flg3, fac1, fac2, fac3, v, sv, st_idx)
    end
end
#--------------------------------------------------------
function interpolate_brick_CUDA!(offset::AbstractArray{T,1},  m1_r::AbstractArray{FT,1},   wb::AbstractArray{FT,1}, 
                                     ξ1::AbstractArray{FT,1},   ξ2::AbstractArray{FT,1},   ξ3::AbstractArray{FT,1}, 
                                     o1::AbstractArray{T,1},    o2::AbstractArray{T,1},    o3::AbstractArray{T,1}, 
                                   flg1::AbstractArray{T,1},  flg2::AbstractArray{T,1},  flg3::AbstractArray{T,1}, 
                                   fac1::AbstractArray{FT,1}, fac2::AbstractArray{FT,1}, fac3::AbstractArray{FT,1}, 
                                      v::AbstractArray{FT,1},   sv::AbstractArray{FT}, st_idx::T) where {T <: Int, FT <: AbstractFloat}



    tj = threadIdx().x; tk = threadIdx().y; # thread ids	
    el = blockIdx().x                       # assigning one element per block 

    qm1 = length(m1_r)
    #--------creating views for shared memory
    shm_FT = @cuDynamicSharedMem(FT, (qm1,qm1+2)) 
    
    vout_jk = view(shm_FT,:,1:qm1)
    wb_sh   = view(shm_FT,:,qm1+1) 
    m1_r_sh = view(shm_FT,:,qm1+2) 
    #------loading shared memory-----------------------------
    if tk==1
        wb_sh[tj]   = wb[tj]
        m1_r_sh[tj] = m1_r[tj] 
    end 
    sync_threads() 
    #-----------------------------------
    off = offset[el]
    #-----------------------------------
    s1 = o1[el]+1; e1 = o1[el+1]; l1 = e1-s1+1;
    s2 = o2[el]+1; e2 = o2[el+1]; l2 = e2-s2+1;  
   	s3 = o3[el]+1; e3 = o3[el+1]; l3 = e3-s3+1;

    if l1 > 0
        for k in 1:l3, j in 1:l2, i in 1:l1 # interpolating point-by-point
			f1 = flg1[s1+i-1]; fac1l = fac1[s1+i-1]; ξ1l = ξ1[s1+i-1]
			f2 = flg2[s2+j-1]; fac2l = fac2[s2+j-1]; ξ2l = ξ2[s2+j-1]
			f3 = flg3[s3+k-1]; fac3l = fac3[s3+k-1]; ξ3l = ξ3[s3+k-1]
            if f1==0 # applying phir
                @inbounds vout_jk[tj,tk] = sv[1 + (tj-1)*qm1 + (tk-1)*qm1*qm1, st_idx, el] *  wb_sh[1] / (ξ1l-m1_r_sh[1])  
             	for ii in 2:qm1
                    @inbounds vout_jk[tj,tk] += sv[ii + (tj-1)*qm1 + (tk-1)*qm1*qm1, st_idx, el] *  wb_sh[ii] / (ξ1l-m1_r_sh[ii]) 
                end
            else
                @inbounds vout_jk[tj,tk] = sv[f1 + (tj-1)*qm1 + (tk-1)*qm1*qm1, st_idx, el]
            end
            if f2==0 # applying phis
                @inbounds vout_jk[tj,tk] *= (wb_sh[tj] / (ξ2l-m1_r_sh[tj]))
            end
            sync_threads() 

	        if tj==1 # reduction
    	        if f2==0
        	        for ij in 2:qm1
            	        @inbounds vout_jk[1,tk]  += vout_jk[ij,tk]
                	end
	            else
    	            if f2 ≠ 1
        	            @inbounds vout_jk[1,tk] = vout_jk[f2,tk]
            	    end
	            end
	
    	        if f3==0 # applying phit
        	        @inbounds vout_jk[1,tk] *= (wb_sh[tk] / (ξ3l-m1_r_sh[tk]))
            	end
	        end
    	    sync_threads()

	        if tj==1 && tk==1 # reduction
    	        if f3==0
        	        for ik in 2:qm1
            	        @inbounds vout_jk[1,1] += vout_jk[1,ik] 
                	end
	            else
    	            if f3 ≠ 1
        	            @inbounds vout_jk[1,1] = vout_jk[1,f3]
            	    end
	            end
    	        @inbounds v[off + i + (j-1)*l1 + (k-1)*l1*l2] = ( vout_jk[1,1] * fac1l * fac2l * fac3l )
        	end
        end

    end

    return nothing
end
#--------------------------------------------------------
"""
    interpolationvector_1pt!(rsrc::Vector{FT}, rdst::FT,
                             wbsrc::Vector{FT}, phi::Vector{FT}) where FT <: AbstractFloat

returns the polynomial interpolation matrix for interpolating between the points
`rsrc` (with associated barycentric weights `wbsrc`) and a single point rdst. This 
function writes the result to a vector since the basis functions are calculated at a single point

Reference:
  Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange Interpolation",
  SIAM Review 46 (2004), pp. 501-517.
  <https://doi.org/10.1137/S0036144502417715>
"""
@inline function interpolationvector_1pt!(rsrc::AbstractVector{FT}, rdst::FT,
                             wbsrc::AbstractVector{FT}, phi::AbstractVector{FT}) where FT <: AbstractFloat
    qm1 = length(rsrc)
    @assert length(phi) == qm1

    for ib in 1:qm1
        if rdst==rsrc[ib]
            phi .= FT(0)
            @inbounds phi[ib] = FT(1) 
            break
        else
            @inbounds phi[ib] = wbsrc[ib] / (rdst-rsrc[ib])
        end
    end
    d = sum(phi)
    phi ./= d 
    return nothing
end
#--------------------------------------------------------
"""
    InterpolationCubedSphere(grid::DiscontinuousSpectralElementGrid, vert_range::AbstractArray{FT}, lat_res::FT, long_res::FT, rad_res::FT) where {FT <: AbstractFloat}

This interpolation structure and the corresponding functions works for a cubed sphere topology. The data is interpolated along a lat/long/rad grid.

# input for the inner constructor
 - `grid` DiscontinousSpectralElementGrid
 - `vert_range` vertex range along the radial coordinate 
 - `lat_res` Resolution of the interpolation grid along the latitude coordinate in radians 
 - `long_res` Resolution of the interpolation grid along the longitude coordinate in radians 
 - `rad_res` Resolution of the interpolation grid along the radial coordinate 
"""
struct InterpolationCubedSphere{T <: Int, FT <: AbstractFloat, TA<:AbstractArray{T,1}, FTA<:AbstractArray{FT,1}}

    realelems::UnitRange{Int64}
    poly_order::T

    lat_min::FT;  long_min::FT;  rad_min::FT; # domain bounds, min
    lat_max::FT;  long_max::FT;  rad_max::FT; # domain bounds, max
    lat_res::FT;  long_res::FT;  rad_res::FT; # respective resolutions for the uniform grid

    n_lat::T; n_long::T; n_rad::T;            # # of lat, long & rad grid locations

    ξ1::FTA #Device array containing ξ1 coordinates of interpolation points within each element 
    ξ2::FTA #Device array containing ξ2 coordinates of interpolation points within each element 
    ξ3::FTA #Device array containing ξ3 coordinates of interpolation points within each element 

    flg::AbstractArray{T,2} # flags when ξ1/ξ2/ξ3 interpolation point matches with a GLL point

    fac::FTA # normalization factor

    radc::FTA  # rad coordinates of interpolation points within each element
    latc::FTA  # lat coordinates of interpolation points within each element
    longc::FTA # long coordinates of interpolation points within each element


    v::FTA      # interpolated variable within each element
    offset::TA  # offsets for each element for v 

    m1_r::FTA   # GLL points
    m1_w::FTA   # GLL weights
    wb::FTA     # Barycentric weights

  #--------------------------------------------------------
    function InterpolationCubedSphere(grid::DiscontinuousSpectralElementGrid, vert_range::AbstractArray{FT}, nhor::Int, lat_res::FT, long_res::FT, rad_res::FT) where {FT <: AbstractFloat}
		DA = CLIMA.array_type()                    # device array
        poly_order = polynomialorder(grid)
        qm1    = poly_order + 1
        toler1 = FT(eps(FT) * vert_range[1] * 2.0) # tolerance for unwarp function
        toler2 = FT(eps(FT) * 4.0)                 # tolerance 
        toler3 = FT(eps(FT) * vert_range[1] * 10.0) # tolerance for Newton-Raphson 

        lat_min,   lat_max = FT(0.0), FT(π)                 # inclination/zeinth angle range
        long_min, long_max = FT(0.0), FT(2*π)  			    # azimuthal angle range
        rad_min,   rad_max = vert_range[1], vert_range[end] # radius range

        realelems = grid.topology.realelems # Element (numbers) on the local processor
        Nel = length(realelems)

        nvert = length(vert_range) - 1              # # of elements in vertical direction
        Nel_glob = nvert * nhor * nhor * 6

        nblck = nhor * nhor * nvert
        Δh = 2 / nhor                               # horizontal grid spacing in unwarped grid

        lat_grd, long_grd, rad_grd = range(lat_min, lat_max, step=lat_res), range(long_min, long_max, step=long_res), range(rad_min, rad_max, step=rad_res) 

        n_lat, n_long, n_rad = Int(length(lat_grd)), Int(length(long_grd)), Int(length(rad_grd))

        uw_grd = zeros(FT, 3, 1)
        #---------------------------------------------- 
        glob_ord = grid.topology.origsendorder # to account for reordering of elements after the partitioning process 

        glob_elem_no = zeros(Int, nvert*length(glob_ord))

        for i in 1:length(glob_ord), j in 1:nvert
            glob_elem_no[j + (i-1)*nvert] = (glob_ord[i] - 1)*nvert + j 
        end

        ξ1, ξ2, ξ3        = map( i -> zeros(FT,i), zeros(Int,Nel)), map( i -> zeros(FT,i), zeros(Int,Nel)), map( i -> zeros(FT,i), zeros(Int,Nel))
        radc, latc, longc = map( i -> zeros(FT,i), zeros(Int,Nel)), map( i -> zeros(FT,i), zeros(Int,Nel)), map( i -> zeros(FT,i), zeros(Int,Nel))

        offset = Vector{Int}(undef,Nel+1)
        offset[1] = 0


        for k in 1:n_long, j in 1:n_lat, i in 1:n_rad 
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
                vert_range[1] - rad < toler1 ? l_nrm = 1 :  error("fatal error, rad lower than inner radius: ", vert_range[1] - rad," $x1_grd /// $x2_grd //// $x3_grd" )
            elseif rad ≥ vert_range[end] # accounting for minor rounding errors from unwarp function at boundaries 
                rad - vert_range[end] < toler1 ? l_nrm = nvert : error("fatal error, rad greater than outer radius")
            else                         # normal scenario
                l_nrm = findfirst( X -> X .- rad .> 0.0, vert_range ) - 1 # identify stack bin 
            end
            #--------------------------------
            if     abs(x1_uw2_grd + 1) < toler2 # face 1 (x1 == -1 plane)
	    	    l2 = min(div(x2_uw2_grd + 1, Δh) + 1, nhor)
    		    l3 = min(div(x3_uw2_grd + 1, Δh) + 1, nhor)
                el_glob = Int(l_nrm + (nhor-l2)*nvert + (l3-1)*nvert*nhor)
            elseif abs(x2_uw2_grd + 1) < toler2 # face 2 (x2 == -1 plane)
		        l1 = min(div(x1_uw2_grd + 1, Δh) + 1, nhor)
    		    l3 = min(div(x3_uw2_grd + 1, Δh) + 1, nhor)
                el_glob = Int(l_nrm + (l1-1)*nvert + (l3-1)*nvert*nhor + nblck*1)
            elseif abs(x1_uw2_grd - 1) < toler2 # face 3 (x1 == +1 plane)
		        l2 = min(div(x2_uw2_grd + 1, Δh) + 1, nhor)
    		    l3 = min(div(x3_uw2_grd + 1, Δh) + 1, nhor)
                el_glob = Int(l_nrm + (l2-1)*nvert + (l3-1)*nvert*nhor + nblck*2 )
            elseif abs(x3_uw2_grd - 1) < toler2 # face 4 (x3 == +1 plane)
      		    l1 = min(div(x1_uw2_grd + 1, Δh) + 1, nhor)
	    	    l2 = min(div(x2_uw2_grd + 1, Δh) + 1, nhor)
                el_glob = Int(l_nrm + (l1-1)*nvert + (l2-1)*nvert*nhor + nblck*3)
            elseif abs(x2_uw2_grd - 1) < toler2 # face 5 (x2 == +1 plane)
	            l1 = min(div(x1_uw2_grd + 1, Δh) + 1, nhor)
		        l3 = min(div(x3_uw2_grd + 1, Δh) + 1, nhor)
                el_glob = Int(l_nrm + (l1-1)*nvert + (nhor-l3)*nvert*nhor + nblck*4 )
            elseif abs(x3_uw2_grd + 1) < toler2 # face 6 (x3 == -1 plane)
	    	    l1 = min(div(x1_uw2_grd + 1, Δh) + 1, nhor)
		        l2 = min(div(x2_uw2_grd + 1, Δh) + 1, nhor)
                el_glob = Int(l_nrm + (l1-1)*nvert + (nhor-l2)*nvert*nhor + nblck*5)
            else
                error("error: unwrapped grid does not lie on any of the 6 faces")
            end
            #--------------------------------
            el_loc = findfirst(X -> X-el_glob == 0, glob_elem_no)
            if ( el_loc ≠ nothing ) # computing inner coordinates for local elements
                ξ = invert_trilear_mapping_hex(grid.topology.elemtocoord[1,:,el_loc], 
                                               grid.topology.elemtocoord[2,:,el_loc], 
                                               grid.topology.elemtocoord[3,:,el_loc], uw_grd, toler3)
                push!(ξ1[el_loc],ξ[1]) 
                push!(ξ2[el_loc],ξ[2]) 
                push!(ξ3[el_loc],ξ[3]) 
                push!(radc[el_loc],  rad_grd[i])
                push!(latc[el_loc],  lat_grd[j]) 
                push!(longc[el_loc],long_grd[k])
            end
            #--------------------------------
        end

        [ offset[el+1] = offset[el] + length(radc[el]) for el in 1:Nel]

        Np = offset[Nel+1]
 
        v = Vector{FT}(undef,offset[Nel+1]) # Allocating storage for interpolation variable

        ξ1_d = Array{FT}(undef,Np); ξ2_d = Array{FT}(undef,Np); ξ3_d = Array{FT}(undef,Np)

        flg_d = zeros(Int,3,Np); fac_d  = ones(FT,Np)

        rad_d  = Array{FT}(undef,Np); lat_d  = Array{FT}(undef,Np); long_d = Array{FT}(undef,Np)

        m1_r, m1_w = GaussQuadrature.legendre(FT,qm1,GaussQuadrature.both)
        wb = Elements.baryweights(m1_r)


        for i in 1:Nel
            ctr = 1
            for j in offset[i]+1:offset[i+1]
                ξ1_d[j]    = ξ1[i][ctr]
                ξ2_d[j]    = ξ2[i][ctr]
                ξ3_d[j]    = ξ3[i][ctr]
                rad_d[j]  = radc[i][ctr]
                lat_d[j]  = latc[i][ctr]
                long_d[j] = longc[i][ctr]
                #-----setting up interpolation 
                l1 = findfirst(X -> abs.(X.-ξ1_d[j]) < toler2, m1_r); fac1 = FT(0)
                l2 = findfirst(X -> abs.(X.-ξ2_d[j]) < toler2, m1_r); fac2 = FT(0)
                l3 = findfirst(X -> abs.(X.-ξ3_d[j]) < toler2, m1_r); fac3 = FT(0)

                l1 == nothing ? [ fac1 +=  wb[ib] / (ξ1_d[j]-m1_r[ib]) for ib in 1:qm1 ] : (flg_d[1,j] = l1; fac1 = FT(1);)
                l2 == nothing ? [ fac2 +=  wb[ib] / (ξ2_d[j]-m1_r[ib]) for ib in 1:qm1 ] : (flg_d[2,j] = l2; fac2 = FT(1);)
                l3 == nothing ? [ fac3 +=  wb[ib] / (ξ3_d[j]-m1_r[ib]) for ib in 1:qm1 ] : (flg_d[3,j] = l3; fac3 = FT(1);)

                fac_d[j] = FT(1) / (fac1 * fac2 * fac3)

                #-----------------------------
                ctr += 1
            end 
        end


        ξ1_d = DA(ξ1_d); ξ2_d = DA(ξ2_d); ξ3_d = DA(ξ3_d)

        flg_d = DA(flg_d); fac_d = DA(fac_d)

        rad_d = DA(rad_d); lat_d = DA(lat_d); long_d = DA(long_d)

        m1_r = DA(m1_r); m1_w = DA(m1_w); wb = DA(wb);
        
        offset = DA(offset); v = DA(v);
println("DA = ", DA)
println("typeof(v) = ", typeof(v))
println("typeof(offset) = ", typeof(offset))

        return new{Int, FT, AbstractArray{Int,1}, AbstractArray{FT,1}}(realelems, poly_order, lat_min, long_min, rad_min, lat_max, long_max, rad_max, lat_res, long_res, rad_res, 
                    n_lat, n_long, n_rad, ξ1_d, ξ2_d, ξ3_d, flg_d, fac_d, rad_d, lat_d, long_d, v, offset, m1_r, m1_w, wb)

    #-----------------------------------------------------------------------------------
    end # Inner constructor function InterpolationCubedSphere
#-----------------------------------------------------------------------------------
end # structure InterpolationCubedSphere
#--------------------------------------------------------
"""
    invert_trilear_mapping_hex(X1::Array{FT}, X2::Array{FT}, X3::Array{FT}, x::Array{FT}, tol::FT) where FT <: AbstractFloat 

This function computes (ξ1,ξ2,ξ3) given (x1,x2,x3) and the (8) vertex coordinates of a Hexahedron. Newton-Raphson method is used
# input
 - `X1` X1 coordinates of the (8) vertices of the hexahedron
 - `X2` X2 coordinates of the (8) vertices of the hexahedron
 - `X3` X3 coordinates of the (8) vertices of the hexahedron
 - `x` (x1,x2,x3) coordinates of the 
"""
function invert_trilear_mapping_hex(X1::Array{FT}, X2::Array{FT}, X3::Array{FT}, x::Array{FT}, tol::FT) where FT <: AbstractFloat 
    max_it = 10        # maximum # of iterations
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
            error("invert_trilinear_mapping_hex: Newton-Raphson not converging to desired tolerance after max_it = ", max_it," iterations; err = ", err,"; toler = ", tol)
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
    Jac = MMatrix{3,3,FT,9}(undef)
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
function interpolate_cubed_sphere!(offset::AbstractArray{T,1}, m1_r::AbstractArray{FT,1}, wb::AbstractArray{FT,1}, 
                                   ξ1::AbstractArray{FT,1}, ξ2::AbstractArray{FT,1}, ξ3::AbstractArray{FT,1}, 
                                   flg::AbstractArray{T,2}, fac::AbstractArray{FT,1},
                                   v::AbstractArray{FT,1}, sv::AbstractArray{FT}, st_no::T) where {T <: Integer, FT <: AbstractFloat}

    qm1 = length(m1_r)
    Nel = length(offset) - 1

	device = typeof(sv) <: Array ? CPU() : CUDA()

    if device==CPU()
println("calling CPU")
        #------------------------------------------------------------------------------------------
        Nel = length(offset) - 1

        vout    = FT(0)
        vout_ii = FT(0)
        vout_ij = FT(0)

        for el in 1:Nel #-----for each element elno 
            np  = offset[el+1] - offset[el]
            lag = @view sv[:,st_no,el]
            off = offset[el]

            for i in 1:np # interpolating point-by-point
                ξ1l = ξ1[off+i]; ξ2l = ξ2[off+i]; ξ3l = ξ3[off+i]

                f1 = flg[1,off+i]; f2 = flg[2,off+i]; f3 = flg[3,off+i]
                fc = fac[off+i]

                vout = 0.0
                f3 == 0 ? (ikloop = 1:qm1) : (ikloop = f3:f3)

                for ik in ikloop
                    #--------------------------------------------
                    vout_ij = 0.0
                    f2 == 0 ? (ijloop = 1:qm1) : (ijloop = f2:f2)
                    for ij in ijloop #1:qm1 
                        #----------------------------------------------------------------------------
                        vout_ii = 0.0 

                        if f1 == 0
                            for ii in 1:qm1
                                @inbounds vout_ii += lag[ii + (ij-1)*qm1 + (ik-1)*qm1*qm1] *  wb[ii] / (ξ1l-m1_r[ii])#phir[ii]
                            end
                        else
                            @inbounds vout_ii = lag[f1 + (ij-1)*qm1 + (ik-1)*qm1*qm1]
                        end
                        if f2==0
                            @inbounds vout_ij += vout_ii * wb[ij] / (ξ2l-m1_r[ij])#phis[ij]
                        else
                            @inbounds vout_ij = vout_ii
                        end
                        #----------------------------------------------------------------------------
                    end
                    if f3==0
                        @inbounds vout += vout_ij * wb[ik] / (ξ3l-m1_r[ik])#phit[ik]
                    else
                        @inbounds vout = vout_ij
                    end
                    #--------------------------------------------
                end
                @inbounds v[off + i] = vout * fc

            end
        end
        #------------------------------------------------------------------------------------------
    else
        #------------------------------------------------------------------------------------------
println("calling CUDA")
        @cuda threads=(qm1,qm1) blocks=Nel shmem=qm1*(qm1+2)*sizeof(FT) interpolate_cubed_sphere_CUDA!(offset, m1_r, wb, ξ1, ξ2, ξ3, flg, fac, v, sv, st_no)
        #------------------------------------------------------------------------------------------
    end
    return nothing
end
#--------------------------------------------------------
function interpolate_cubed_sphere_CUDA!(offset::AbstractArray{T,1}, m1_r::AbstractArray{FT,1}, wb::AbstractArray{FT,1}, 
                                        ξ1::AbstractArray{FT,1}, ξ2::AbstractArray{FT,1}, ξ3::AbstractArray{FT,1}, 
                                        flg::AbstractArray{T,2}, fac::AbstractArray{FT,1},
                                        v::AbstractArray{FT,1}, sv::AbstractArray{FT}, st_no::T) where {T <: Integer, FT <: AbstractFloat}

    tj = threadIdx().x; tk = threadIdx().y; # thread ids	
    el = blockIdx().x                       # assigning one element per block 


    qm1 = length(m1_r)
    #--------creating views for shared memory
    shm_FT = @cuDynamicSharedMem(FT, (qm1,qm1+2)) 
    
    vout_jk = view(shm_FT,:,1:qm1)
    wb_sh   = view(shm_FT,:,qm1+1) 
    m1_r_sh = view(shm_FT,:,qm1+2) 
    #------loading shared memory-----------------------------
    if tk==1
        wb_sh[tj]   = wb[tj]
        m1_r_sh[tj] = m1_r[tj] 
    end 
    sync_threads() 
    #-----------------------------------


    #-----------------------------------
    np  = offset[el+1] - offset[el]
    off = offset[el]

    for i in 1:np # interpolating point-by-point
        #-----------------------------------
        ξ1l = ξ1[off+i]; ξ2l = ξ2[off+i]; ξ3l = ξ3[off+i]

        f1 = flg[1,off+i]; f2 = flg[2,off+i]; f3 = flg[3,off+i]
        fc = fac[off+i]


        if f1==0 # applying phir
            @inbounds vout_jk[tj,tk] = sv[1 + (tj-1)*qm1 + (tk-1)*qm1*qm1, st_no, el] *  wb_sh[1] / (ξ1l-m1_r_sh[1])  
           	for ii in 2:qm1
                @inbounds vout_jk[tj,tk] += sv[ii + (tj-1)*qm1 + (tk-1)*qm1*qm1, st_no, el] *  wb_sh[ii] / (ξ1l-m1_r_sh[ii]) 
            end
        else
            @inbounds vout_jk[tj,tk] = sv[f1 + (tj-1)*qm1 + (tk-1)*qm1*qm1, st_no, el]
        end

        if f2==0 # applying phis
            @inbounds vout_jk[tj,tk] *= (wb_sh[tj] / (ξ2l-m1_r_sh[tj]))
        end
        sync_threads() 
            
        if tj==1 # reduction
            if f2==0
                for ij in 2:qm1
                    @inbounds vout_jk[1,tk]  += vout_jk[ij,tk]
                end
            else
                if f2 ≠ 1
                    @inbounds vout_jk[1,tk] = vout_jk[f2,tk]
                end
            end

            if f3==0 # applying phit
                @inbounds vout_jk[1,tk] *= (wb_sh[tk] / (ξ3l-m1_r_sh[tk]))
            end
        end
        sync_threads()

        if tj==1 && tk==1 # reduction
            if f3==0
                for ik in 2:qm1
                    @inbounds vout_jk[1,1] += vout_jk[1,ik] 
                end
            else
                if f3 ≠ 1
                    @inbounds vout_jk[1,1] = vout_jk[1,f3]
                end
            end
            @inbounds v[off + i] = vout_jk[1,1] * fc
        end


            #-----------------------------------
    end
    #-----------------------------------

    return nothing
end
#--------------------------------------------------------
end # module interploation
