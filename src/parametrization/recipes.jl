# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(lpars::LocalParameters, D::SpatialData)
  @assert ndims(lpars) == 2 "plot only available for 2D local parameters"
  quats = lpars.rotation
  u = [quat2dcm(quats[x])[1,1] for x in 1:length(quats)]
  v = [quat2dcm(quats[x])[1,2] for x in 1:length(quats)]
  x = [coordinates(centroid(D,x))[1] for x in 1:nelements(D)]
  y = [coordinates(centroid(D,x))[2] for x in 1:nelements(D)]
  c = lpars.magnitude[2,:]

  seriestype --> :quiver
  quiver --> (u,v)
  line_z --> repeat(c, inner=4)
  seriescolor --> :redsblues

  x,y
end
