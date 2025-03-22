# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# 2d arrows at Plots.jl 

@recipe function f(lpars::LocalAnisotropy, D::SpatialData, axis=:Y)
  @assert ndims(lpars) == 2 "plot only available for 2D local anisotropies"
  quats = lpars.rotation
  u = [quat_to_dcm(quats[x])[1, 1] for x in 1:length(quats)]
  v = [quat_to_dcm(quats[x])[1, 2] for x in 1:length(quats)]
  x = [to(centro(D, x))[1] for x in 1:nvals(D)]
  y = [to(centro(D, x))[2] for x in 1:nvals(D)]
  c = lpars.magnitude[iaxis(axis), :]

  seriestype --> :quiver
  quiver --> (u, v)
  line_z --> repeat(c, inner=4)
  if maximum(c) <= 1
    seriescolor --> :redsblues
  else
    seriescolor --> :bluesreds
  end

  x, y
end
