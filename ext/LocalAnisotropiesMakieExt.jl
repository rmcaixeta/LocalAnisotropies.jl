# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module LocalAnisotropiesMakieExt

using LocalAnisotropies
import Makie

Makie.@recipe(LA, lpars, geom) do scene
    return Makie.Attributes(;
        alpha = 0.6,
		size  = 0.5,
		color = :black,
		colormap = :RdYlBu,
    )
end

Makie.plottype(::LocalAnisotropy, ::LocalAnisotropies.SpatialData) = LA{<:Tuple{LocalAnisotropy, LocalAnisotropies.SpatialData}}

function pre_ellipsoid(lpars, geom, f)
  # pre process
  coords = LocalAnisotropies.coords(geom)
  D3 = size(coords,1) >= 3
  x = coords[1,:]
  y = coords[2,:]
  z = D3 ? coords[3,:] : map(i->0, x)
  s = D3 ? [f * Makie.Vec3f(i[1],i[2],i[3]) for i in eachcol(lpars.magnitude)] : [f * Makie.Vec3f(i[1],i[2],i[2]) for i in eachcol(lpars.magnitude)]
  q = map(q->Makie.Quaternion(q.q1, q.q2, q.q3, q.q0), lpars.rotation)

  x, y, z, s, q
end

function Makie.convert_arguments(P::Type{<:LA}, lpars::LocalAnisotropy, geom::LocalAnisotropies.SpatialData)
	lpars, geom
end

"""
    Makie.plot(localaniso, geometry; alpha=0.6, size=0.5, color=:black)
    Makie.plot!(localaniso, geometry; alpha=0.6, size=0.5, color=:black)

Plot local anisotropies as ellipsoid in Makie.jl environment. Example in the
README file. If `color` is a vector, `colormap` can be defined (default=`:RdYlBu`)
"""
function Makie.plot!(plot::LA)
  # input pars
  lpars = plot[1]
  geom  = plot[2]
  size  = plot[:size]

  out = Makie.@lift pre_ellipsoid($lpars, $geom, $size)
  x = Makie.@lift $out[1]
  y = Makie.@lift $out[2]
  z = Makie.@lift $out[3]
  s = Makie.@lift $out[4]
  q = Makie.@lift $out[5]

  Makie.meshscatter!(plot, x, y, z, markersize=s, rotations=q,
     alpha=plot[:alpha], color=plot[:color], colormap=plot[:colormap])
end

end
