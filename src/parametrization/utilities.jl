# make rescale function more customizable for ratio/ranges uses
# export to vtk scaled ellipses/ellipsoid

"""
    rescale_magnitude(localpars, r1, r2=nothing, r3=nothing;
    clip=[0.05,0.95], magnitude=:ratios)

Rescale magnitude values of `localpars` to desired limits using min-max scaling.
`r1` and `r2` set the limits for the ratios between (range2/range1) and
(range3/range1) respectively. `clip` define the quantiles to which min and max
limits are defined. The default is set to [0.05,0.95] to avoid outliers influence.

## Example

In the 2-D example, the max local anisotropy is 5x and the minimum is 1x
(isotropic). In the 3-D example, the max local anisotropy is 2x and the minimum is 1x

```julia
example2d = rescale_magnitude(localpars2d, (0.2,1.0))
example3d = rescale_magnitude(localpars3d, (0.5,1.0), (0.1,0.5))
```
"""
function rescale_magnitude(lp::LocalParameters, r1, r2=nothing, r3=nothing;
	                       clip=[0.05,0.95], magnitude=:ratios)
    N = ndims(lp)
    m = lp.magnitude
    bounds = [r1,r2]

    for i in 2:N
		bounds[i-1] == nothing && continue
        mi = m[i,:]
        qx = quantile(mi, clip)
        mi[mi .< qx[1]] .= qx[1]
        mi[mi .> qx[2]] .= qx[2]
        r = (mi .- qx[1]) ./ (qx[2]-qx[1])
        typeof(bounds[i-1]) <: Number && (bounds[i-1]=(bounds[i-1],bounds[i-1]))
        b1, b2 = bounds[i-1]
        m[i,:] .= (b1 .+ (b2-b1) .* r)
    end

    LocalParameters(lp.rotation,m)
end

"""
    localpars2vtk(vtkfile, coords, localpars; dir=:all, magnitude=:ratios)

Export local parameters `localpars` at `coords` to VTK format. It export all
directions by default, unless `dir` is set to `:X`, `:Y` or `:Z`. The magnitude
can be expressed as `:ratios` (default) or `:ranges`. In the case of `:ratios`,
`ratio1=range2/range1` and `ratio2=range3/range1`. Instead of an array of
coordinates, `coords` can also be just the georeferenced spatial object. The
output `.vtu` file can be loaded in Paraview or similars and processed as Glyphs
to plot the directions.
"""
function localpars2vtk(vtkfile, coords::AbstractArray, lpars; dir=:all, magnitude=:ratios)
	@assert dir in [:all, :X, :Y, :Z] "`dir` must be :all, :X, :Y or :Z"
	@assert magnitude in [:ratios, :ranges] "`magnitude` must be :ratios or :ranges"

	cod = Dict(:X=>1,:Y=>2,:Z=>3)
	dim = size(coords,1)
	n   = size(coords,2)
	axs = dir == :all ? (1:dim) : [cod[dir]]
	xyz = dim == 2 ? vcat(coords, zeros(Float64,1,n)) : coords
	vtx = [MeshCell( VTKCellTypes.VTK_VERTEX, [i]) for i in 1:n]
	ijk = [zeros(Float64, 3, n, 1, 1) for i in axs]
	for x in 1:n
		dcm = quat_to_dcm(lpars.rotation[x])
        f = dcm[3,3] < 0 ? -1 : 1
		for dir in axs
			ijk[dir][1, x, 1, 1] = f*dcm[dir,1]
			ijk[dir][2, x, 1, 1] = f*dcm[dir,2]
			ijk[dir][3, x, 1, 1] = f*dcm[dir,3]
		end
	end
	outfiles = vtk_grid(vtkfile,xyz,(vtx)) do vtk
		for dir in axs
			vtk["dir$dir"] = ijk[dir]
		end
		for dir in 1:dim
			if magnitude == :ranges
				vtk["range$dir"] = lpars.magnitude[dir,:]
			elseif dir > 1
				vtk["ratio$(dir-1)"] = lpars.magnitude[dir,:] ./ lpars.magnitude[1,:]
			end
		end
	end
end

function localpars2vtk(vtkfile, D::SpatialData, lpars; kwargs...)
	coords = reduce(hcat, [coordinates(centroid(D,x)) for x in 1:nelements(D)])
	localpars2vtk(vtkfile, coords, lpars; kwargs...)
end

toqmat(lp) = [qmat(rotation(lp,i),magnitude(lp,i)) for i in 1:nvals(lp)]
