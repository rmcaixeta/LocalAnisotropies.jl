# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    localvariography(sdata, lpars, var₁, var₂=var₁;
    axis=:X, tol=1e-6, maxratio1=Inf, maxratio2=Inf, [optional parameters])

Computes an unconventional directional (cross-)variogram for the variables `var₁`
and `var₂` of spatial data `sdata`. This uses local anisotropies `lpars` to create
pseudo-partitions along local `axis` direction (`axis` ∈ [:X, :Y, :Z]) and the
EmpiricalVariogram is calculated with these data. The directional tolerance
bandwidth `tol` can be passed. The `maxratio1` and `maxratio2` partially ignore
data if magnitude ratios are above these values (does not ignore anything by
default). Other optional parameters are defined at `GeoStats.EmpiricalVariogram`,
such as `nlags`, `maxlag`, `distance` and `algo`.

Similar (but not equal) to https://github.com/rmcaixeta/Local_variography
"""
function localvariography(obj::SpatialData, lpars::LocalAnisotropy, var₁, var₂=var₁;
                          axis=:X, tol=1e-6, maxratio1=Inf, maxratio2=Inf, kwargs...)
  p = pseudolocalpartition(obj, lpars, axis, tol, maxratio1, maxratio2)
  EmpiricalVariogram(p, var₁, var₂; kwargs...)
end

function pseudolocalpartition(obj, lpars, axis, tol, maxratio1, maxratio2)
  subs = []
  dims = ndims(lpars)
  n = nvals(obj)
  for i in 1:n
    maxratio1 != Inf && ratio1(lpars,i) > maxratio1 && continue
    maxratio2 != Inf && ratio2(lpars,i) > maxratio2 && continue
    x = coordinates(centroid(obj, i))
    v = rotmat(lpars, i)[iaxis(axis),1:dims]
    p = DirectionPartition(Tuple(v), tol=tol)
    s = Int[i]
    for j in 1:n
      i == j && continue
      y = coordinates(centroid(obj, j))
      p(x, y) && push!(s, j)
    end
    push!(subs, s)
  end

  Partition(obj, subs)
end
