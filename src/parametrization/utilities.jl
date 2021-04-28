# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    rescale_magnitude(localaniso, r1, r2=nothing, r3=nothing;
    clip=[0.05,0.95], magnitude=:ratios)
	rescale_magnitude(localaniso; r1=nothing, r2=nothing, r3=nothing
    clip=[0.05,0.95], magnitude=:ratios)

Rescale magnitude values of `localaniso` to desired limits using min-max scaling.
`r1` and `r2` set the limits for the ratios between (range2/range1) and
(range3/range1) respectively. `clip` define the quantiles to which min and max
limits are defined. The default is set to [0.05,0.95] to avoid outliers influence.

## Example

In the 2-D example, the max local anisotropy is 5x and the minimum is 1x
(isotropic). In the 3-D example, the max local anisotropy is 2x and the minimum is 1x

```julia
example2d = rescale_magnitude(localaniso2d, (0.2,1.0))
example3d = rescale_magnitude(localaniso3d, r1=(0.5,1.0), r2=(0.1,0.5))
```
"""
function rescale_magnitude!(lp::LocalAnisotropy, r1, r2=nothing, r3=nothing;
	                       clip=[0.05,0.95], magnitude=:ratios)
    N = ndims(lp)
    m = lp.magnitude
	rvals = [r1,r2,r3]

	if magnitude == :ratios
		ref = m[1,:]
		fixaxis!(lp, :X)
		rvals = [nothing,r1,r2]
	end

    for i in 1:N
		rvals[i] == nothing && continue
        mi = m[i,:]
        qx = quantile(mi, clip)
        mi[mi .< qx[1]] .= qx[1]
        mi[mi .> qx[2]] .= qx[2]
        mi = (mi .- qx[1]) ./ (qx[2]-qx[1])
        rvals[i] isa Number && (rvals[i]=(rvals[i],rvals[i]))
        b1, b2 = rvals[i]
        m[i,:] .= (b1 .+ (b2-b1) .* mi)
		magnitude == :ratios && (m[i,:] .*= ref)
    end
	lp
end

function rescale_magnitude!(lp::LocalAnisotropy; r1=nothing, r2=nothing,
	                       r3=nothing, kwargs...)
    rescale_magnitude!(lp, r1, r2, r3; kwargs...)
end

function rescale_magnitude(lpo::LocalAnisotropy, r1, r2=nothing, r3=nothing;
	                       clip=[0.05,0.95], magnitude=:ratios)
    lp = deepcopy(lpo)
	rescale_magnitude!(lp, r1, r2, r3; clip=clip, magnitude=magnitude)
end

function rescale_magnitude(lpo::LocalAnisotropy; r1=nothing, r2=nothing,
	                       r3=nothing, kwargs...)
    lp = deepcopy(lpo)
	rescale_magnitude!(lp, r1, r2, r3; kwargs...)
end

function fixaxis(lpar::LocalAnisotropy, ax::Symbol)
  lp = deepcopy(lpar)
  fixaxis!(lp, ax)
end

function fixaxis!(lpar::LocalAnisotropy, ax::Symbol)
  # rescale magnitude according to reference axis
  m = lpar.magnitude
  ref = m[iaxis(ax),:]
  for i in 1:size(m, 1)
    m[i,:] ./= ref
  end
  lpar
end

"""
    localaniso2vtk(vtkfile, coords, localaniso; dir=:ellips, magnitude=:ratios)

Export local anisotropies `localaniso` at `coords` to VTK format. It export local
ellipsoids by default. If `dir` is set to `:X`, `:Y`, `:Z` or `:XYZ`, local
vectors/arrows are exported on these axes. The magnitude can be expressed as
`:ratios` (default) or `:ranges`. In the case of `:ratios`, `ratio1=range2/range1`
and `ratio2=range3/range1`. Instead of an array of coordinates, `coords` can
also be just the georeferenced spatial object. The output `.vtu` file can be
loaded in Paraview or similars and processed as Glyphs to plot the directions.
If `dir = :ellips`, TensorGlyphs must be used for proper visualization.
"""
function localaniso2vtk(vtkfile, coords::AbstractArray, lpars; dir=:ellips, magnitude=:ratios)
	@assert dir in [:ellips, :XYZ, :X, :Y, :Z] "`dir` must be :ellips, :XYZ :X, :Y or :Z"
	@assert magnitude in [:ratios, :ranges] "`magnitude` must be :ratios or :ranges"

	elp = (dir == :ellips)
	dim = size(coords,1)
	n   = size(coords,2)
	axs = dir in [:XYZ, :ellips] ? (1:dim) : [iaxis(dir)]
	xyz = dim == 2 ? vcat(coords, zeros(Float64,1,n)) : coords
	vtx = [MeshCell( VTKCellTypes.VTK_VERTEX, [i]) for i in 1:n]
	ijk = elp ? zeros(Float64, 9, n, 1, 1) : [zeros(Float64, 3, n, 1, 1) for i in axs]
	zx  = (dim == 2 && elp) ? minimum(lpars.magnitude)/10 : 0.0

	# put data into vtk structure format
	for x in 1:n
		if elp
			mag = lpars.magnitude[:,x]
			dim == 2 && push!(mag, zx)
			Λ = Diagonal(SVector{3}(mag))
			P = rotmat(lpars, x)' * Λ
			for v in 1:length(P)
				ijk[v, x, 1, 1] = P[v]
			end
		else
			dcm = rotmat(lpars, x)
	        f = dcm[3,3] < 0 ? -1 : 1
			for d in axs
				ijk[d][1, x, 1, 1] = f*dcm[d,1]
				ijk[d][2, x, 1, 1] = f*dcm[d,2]
				ijk[d][3, x, 1, 1] = f*dcm[d,3]
			end
		end
	end

	# write vtk file
	outfiles = vtk_grid(vtkfile,xyz,(vtx)) do vtk
		if elp
			vtk["ellips"] = ijk
		else
			for d in axs
				vtk["dir$d"] = ijk[d]
			end
		end
		for d in 1:dim
			if magnitude == :ranges
				vtk["range$d"] = lpars.magnitude[d,:]
			elseif d > 1
				vtk["ratio$(d-1)"] = lpars.magnitude[d,:] ./ lpars.magnitude[1,:]
			end
		end
	end
end

function localaniso2vtk(vtkfile, D::SpatialData, lpars; kwargs...)
	coords = reduce(hcat, [coordinates(centroid(D,x)) for x in 1:nelements(D)])
	localaniso2vtk(vtkfile, coords, lpars; kwargs...)
end

toqmat(lp) = [qmat(rotation(lp,i),magnitude(lp,i)) for i in 1:nvals(lp)]
