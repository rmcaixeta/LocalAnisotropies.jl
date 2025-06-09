# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    rescale_magnitude(localaniso; r1, r2, r3, clip=[0.05,0.95])
    rescale_magnitude!(localaniso; r1, r2, r3, clip=[0.05,0.95])

Rescale magnitude values of `localaniso` to desired limits using min-max scaling.
`r1`, `r2` and `r3` refer to the ranges along each main axis. `clip` defines
the quantiles to which min and max limits are defined. The default is set to
[0.05,0.95] to avoid outliers influence.

## Example

In the example below, a 3-D local anisotropy object is rescaled. `r1` has all the
values set to 1.0. `r2` and `r3` are rescaled between (0.2, 1.0) and (0.1, 0.5)
respectively. That means that if we associate it to a SphericalVariogram of range=10,
range1 will be equal to 10 at every place. range2 will vary between 2 and 10.
And range3 will vary between 1 and 5.

```julia
example = rescale_magnitude(localaniso3d, r1=1, r2=(0.2,1.0), r3=(0.1,0.5))
```
"""
function rescale_magnitude!(lp::LocalAnisotropy; r1=nothing, r2=nothing, r3=nothing, clip=[0.05, 0.95])
  N = ndims(lp)
  m = lp.magnitude
  rvals = (r1, r2, r3)

  for i in 1:N
    rvals[i] == nothing && continue

    mi = m[i, :]
    qx = quantile(mi, clip)

    if rvals[i] isa Number
      m[i, :] .= rvals[i]
    elseif qx[1] == qx[2]
      throw(ArgumentError("clipped values are unique; inform a single value to r$i"))
    else
      mi[mi .< qx[1]] .= qx[1]
      mi[mi .> qx[2]] .= qx[2]
      mi = (mi .- qx[1]) ./ (qx[2] - qx[1])
      b1, b2 = rvals[i] isa Number ? (rvals[i], rvals[i]) : rvals[i]
      m[i, :] .= (b1 .+ (b2 - b1) .* mi)
    end
  end
  lp
end

function rescale_magnitude(lpo::LocalAnisotropy; kwargs...)
  lp = deepcopy(lpo)
  rescale_magnitude!(lp; kwargs...)
end

"""
    reference_magnitude(localaniso, axis)
    reference_magnitude!(localaniso, axis)

Rescale the magnitudes of `localaniso` by setting `axis` as a unit reference.
`axis` can be :X, :Y or :Z.

## Example

For example, consider a point with magnitude equal to [1.0, 0.2]. After scaling
it with a ExponentialVariogram of range=10m, the range in X would be 10m and in
Y would be 2m at this point. If we first call `reference_magnitude!(localaniso, :Y)`,
the magnitude will turn to [5.0, 1.0], and after scaling to the same previous
variogram, the range in X would be 50m and in Y would be 10m.
"""
function reference_magnitude!(lpar::LocalAnisotropy, ax::Symbol)
  # rescale magnitude according to reference axis
  m = lpar.magnitude
  ref = m[iaxis(ax), :]
  for i in 1:size(m, 1)
    m[i, :] ./= ref
  end
  lpar
end

function reference_magnitude(lpar::LocalAnisotropy, ax::Symbol)
  lp = deepcopy(lpar)
  reference_magnitude!(lp, ax)
end

"""
    adjust_rake(localaniso, azimuths)
    adjust_rake!(localaniso, azimuths)

Adjust the main direction :X along the azimuth informed.
This will change the directions of :X and :Y, and keep the :Z unchanged
Useful when the ellipsoid has the correct :Z (e.g. along body thickness), but the
main direction should align with an specific azimuth onto the main plane.

## Example

Considering three ellipsoids in the `localaniso` variable

```julia
azimuths = [25.0, 30.0, 40.0]
newpars = adjust_rake(localaniso, azimuths)
```
"""
function adjust_rake!(lpar::LocalAnisotropy, az::AbstractVector)
  rot = map(1:nvals(lpar)) do i
    dcm = quat_to_dcm(lpar.rotation[i])
    p1 = Plane(Point(0, 0, 0), Vec(dcm[1, :]), Vec(dcm[2, :]))
    n2 = Vec(sind(az[i] + 90), cosd(az[i] + 90), 0)
    p2 = Plane(Point(0, 0, 0), n2)
    icross = ustrip.(to(intersection(p1, p2).geom.b))
    vectors_to_quaternion(icross, dcm[3, :]; dirs=(1, 3))
  end
  lpar.rotation .= rot
  lpar
end

function adjust_rake(lpar::LocalAnisotropy, az::AbstractVector)
  lp = deepcopy(lpar)
  adjust_rake!(lp, az)
end

adjust_rake!(lpar::LocalAnisotropy, az::Number) = adjust_rake!(lpar, [az for i in 1:nvals(lpar)])
adjust_rake(lpar::LocalAnisotropy, az::Number) = adjust_rake(lpar, [az for i in 1:nvals(lpar)])

"""
    to_vtk(vtkfile, coords, localaniso; dir=:ellips, magnitude=:ranges)

Export local anisotropies `localaniso` at `coords` to VTK format. It export local
ellipsoids by default. If `dir` is set to `:X`, `:Y`, `:Z` or `:XYZ`, local
vectors/arrows are exported on these axes. The magnitude can be expressed as
`:ranges` (default) or `:ratios`. In the case of `:ratios`, `ratio1=range2/range1`
and `ratio2=range3/range1`. Instead of an array of coordinates, `coords` can
also be just the georeferenced spatial object. The output `.vtu` file can be
loaded in Paraview or similars and processed as Glyphs to plot the directions.
If `dir = :ellips`, TensorGlyphs must be used for proper visualization (need to
disable extract eigenvalues).
"""
function to_vtk(vtkfile, coords::AbstractArray, lpars::LocalAnisotropy; dir=:ellips, magnitude=:ranges)
  @assert dir in [:ellips, :XYZ, :X, :Y, :Z] "`dir` must be :ellips, :XYZ :X, :Y or :Z"
  @assert magnitude in [:ratios, :ranges] "`magnitude` must be :ratios or :ranges"

  elp = (dir == :ellips)
  dim = size(coords, 1)
  n = size(coords, 2)
  axs = dir in [:XYZ, :ellips] ? (1:dim) : [iaxis(dir)]
  xyz = dim == 2 ? vcat(coords, zeros(Float64, 1, n)) : coords
  vtx = [MeshCell(VTKCellTypes.VTK_VERTEX, [i]) for i in 1:n]
  ijk = elp ? zeros(Float64, 9, n, 1, 1) : [zeros(Float64, 3, n, 1, 1) for i in axs]
  zx = (dim == 2 && elp) ? minimum(lpars.magnitude) / 10 : 0.0

  # put data into vtk structure format
  for x in 1:n
    if elp
      mag = lpars.magnitude[:, x]
      dim == 2 && push!(mag, zx)
      Λ = Diagonal(SVector{3}(mag))
      P = rotmat(lpars, x)' * Λ
      for v in 1:length(P)
        ijk[v, x, 1, 1] = P[v]
      end
    else
      dcm = rotmat(lpars, x)
      f = dcm[3, 3] < 0 ? -1 : 1
      for d in axs
        ijk[d][1, x, 1, 1] = f * dcm[d, 1]
        ijk[d][2, x, 1, 1] = f * dcm[d, 2]
        ijk[d][3, x, 1, 1] = f * dcm[d, 3]
      end
    end
  end

  # write vtk file
  outfiles = vtk_grid(vtkfile, xyz, (vtx)) do vtk
    if elp
      vtk["ellips"] = ijk
    else
      for d in axs
        vtk["dir$d"] = ijk[d]
      end
    end
    for d in 1:dim
      if magnitude == :ranges
        vtk["range$d"] = lpars.magnitude[d, :]
      elseif d > 1
        vtk["ratio$(d-1)"] = lpars.magnitude[d, :] ./ lpars.magnitude[1, :]
      end
    end
  end
end

function to_vtk(vtkfile, D::SpatialData, lpars::LocalAnisotropy; kwargs...)
  coords = reduce(hcat, [ustrip.(to(centro(D, x))) for x in 1:nvals(D)])
  to_vtk(vtkfile, coords, lpars; kwargs...)
end

function to_vector(lpars::LocalAnisotropy, dir::Int)
  mapreduce(vcat, 1:nvals(lpars)) do x
    to_vector(rotation(lpars, x), dir)
  end
end

function to_vector_angle(lpars::LocalAnisotropy, dir::Int)
  mapreduce(vcat, 1:nvals(lpars)) do x
    to_vector_angle(rotation(lpars, x), dir)
  end
end

function to_plane_angle(lpars::LocalAnisotropy)
  mapreduce(vcat, 1:nvals(lpars)) do x
    to_plane_angle(rotation(lpars, x))
  end
end

function to_plane_dipvector(lpars::LocalAnisotropy)
  mapreduce(vcat, 1:nvals(lpars)) do x
    to_plane_dipvector(rotation(lpars, x))
  end
end

function to_vector(q::Quaternion, dir::Int)
  dcm = quat_to_dcm(q)
  f = dcm[3, 3] < 0 ? -1 : 1
  dvec = f * dcm[dir, :]
  hcat(dvec...)
end

function to_vector_angle(q::Quaternion, dir::Int)
  dcm = quat_to_dcm(q)
  i, j, k = dcm[dir, :]
  plunge = rad2deg(asin(abs(k)))
  azimuth = rad2deg(atan(j, i)) % 360
  hcat(azimuth, plunge)
end

function to_plane_angle(q::Quaternion)
  dcm = quat_to_dcm(q)
  nx, ny, nz = dcm[3, :]
  strike = atan(ny, nx)
  dipdir = (rad2deg(strike + π/2) + 360) % 360
  dip = rad2deg(asin(abs(nz)))
  hcat(dipdir, dip)
end

function to_plane_dipvector(q::Quaternion)
  dcm = quat_to_dcm(q)
  nx, ny, nz = dcm[3, :]
  dip_vector = [nz * nx, nz * ny, -(nx^2 + ny^2)]
  hcat(dip_vector...)
end

function angles_to_vector(azimuth, plunge)
  i = cosd(plunge) * sind(azimuth)
  j = cosd(plunge) * cosd(azimuth)
  k = -sind(plunge)
  [i, j, k]
end

function normal_to_quaternion(normal)
  strike = [-normal[2], normal[1], 0]
  v1 = sum(strike) > 0 ? strike ./ norm(strike) : [1, 0, 0]
  v2 = cross(v1, normal)
  vectors_to_quaternion(v1, v2, normal)
end

function dipvector_to_quaternion(dipdirvec)
    dx, dy, dz = dipdirvec
    hm = (dx^2 + dy^2) ^ 0.5
    n = hm != 0 ? [-dx * dz / hm, -dy * dz / hm, hm] : [1,0,0]
    n ./= norm(n)
    vectors_to_quaternion(dipdirvec, n; dirs=(2, 3))
end

function vectors_to_quaternion(vec1, vec2, vec3)
  m = SMatrix{3,3}(vec1..., vec2..., vec3...)'
  det(m) < 0 && (m = Diagonal(SVector{3}([-1, 1, 1])) * m)
  dcm_to_quat(DCM(m))
end

function vectors_to_quaternion(vec1, vec2; dirs=(1, 2))
  vecs = Vector{Vector{Float64}}(undef, 3)
  vecs[dirs[1]] = vec1
  vecs[dirs[2]] = vec2
  vecs[setdiff(1:3, dirs)[1]] = cross(vec1, vec2)
  vectors_to_quaternion(vecs...)
end

Base.vcat(lpars::LocalAnisotropy...; kwars...) = reduce(vcat, lpars)

function Base.vcat(lpar1::LocalAnisotropy, lpar2::LocalAnisotropy; kind=:union)
  kind != :union && error("vcat of LocalAnisotropy is restricted to kind=:union ")
  rot = vcat(lpar1.rotation, lpar2.rotation)
  mag = hcat(lpar1.magnitude, lpar2.magnitude)
  LocalAnisotropy(rot, mag)
end
