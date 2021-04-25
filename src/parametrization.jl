# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LocalAnisotropy(rotation, magnitude)

Create a `LocalAnisotropy` object. The `rotation` is a vector of
`ReferenceFrameRotations.Quaternion` and the magnitude is an array of size
(number of dimensions, number of elements)
"""
struct LocalAnisotropy
  rotation::AbstractVector{Quaternion}
  magnitude::AbstractArray
end

Ints = Union{Int,AbstractVector{Int}}
rotation(L::LocalAnisotropy) = L.rotation
rotation(L::LocalAnisotropy, i::Int) = L.rotation[i]
rotation(L::LocalAnisotropy, i::AbstractVector{Int}) = view(L.rotation, i)
rotmat(L::LocalAnisotropy, i::Int) = quat_to_dcm(L.rotation[i])
magnitude(L::LocalAnisotropy) = L.magnitude
magnitude(L::LocalAnisotropy, i::Ints) = view(L.magnitude, :, i)
ratio1(L::LocalAnisotropy) = L.magnitude[2,:] ./ L.magnitude[1,:]
ratio1(L::LocalAnisotropy, i::Ints) = L.magnitude[2,i] / L.magnitude[1,i]
ratio2(L::LocalAnisotropy) = L.magnitude[3,:] ./ L.magnitude[1,:]
ratio2(L::LocalAnisotropy, i::Ints) = L.magnitude[3,i] / L.magnitude[1,i]
slice(L::LocalAnisotropy, i::Ints) = LocalAnisotropy(rotation(L,i),magnitude(L,i))
nvals(L::LocalAnisotropy) = length(L.rotation)
ndims(L::LocalAnisotropy) = size(L.magnitude,1)
iaxis(ax) = (X=1, Y=2, Z=3)[ax]

Base.show(io::IO, lp::LocalAnisotropy) = print(io, "LocalAnisotropy")

function Base.show(io::IO, ::MIME"text/plain", lp::LocalAnisotropy)
	print(io,"LocalAnisotropy $(ndims(lp))-D")
end

SpatialData = Union{GeoData,Domain,DataView,DomainView}

struct LocalGeoData
	data::SpatialData
	object::SpatialData
	localaniso::LocalAnisotropy
	refvario::Union{Variogram,Nothing}
	datapars::Union{AbstractVector{Int},Nothing}
    graph::Union{SimpleWeightedGraph,Nothing}
end

"""
    LocalGeoData(object, localaniso)

Association of a spatial data object and its `LocalAnisotropy`.
"""
function LocalGeoData(obj::SpatialData, lpars::LocalAnisotropy)
	LocalGeoData(nothing, obj, lpars, nothing, nothing, nothing)
end

"""
    LocalGeoData(object, localaniso, refvario)

Association of a spatial data object, its `LocalAnisotropy` and its reference variogram.
"""
function LocalGeoData(obj::SpatialData, lpars::LocalAnisotropy, refvario::Variogram)
	LocalGeoData(nothing, obj, lpars, refvario, nothing, nothing)
end

"""
    LocalGeoData(samples, object, localaniso)

Association of samples, a spatial data object and their `LocalAnisotropy`. The
samples `LocalAnisotropy` are passed from the spatial data via nearest neighbors.
"""
function LocalGeoData(hd::SpatialData, obj::SpatialData, lpars::LocalAnisotropy)
	hdids = grid2hd_ids(hd,obj)
	LocalGeoData(hd, obj, lpars, nothing, hdids, nothing)
end

"""
    LocalGeoData(samples, object, localaniso, refvario)

Association of samples, a spatial data object, their `LocalAnisotropy` and the
reference variogram. The samples `LocalAnisotropy` are passed from the spatial data
via nearest neighbors.
"""
function LocalGeoData(hd::SpatialData, obj::SpatialData, lpars::LocalAnisotropy, refvario::Variogram)
	hdids = grid2hd_ids(hd,obj)
	LocalGeoData(hd, obj, lpars, refvario, hdids, nothing)
end

obj(D::LocalGeoData) = D.object
sobj(D::LocalGeoData) = D.data
ndims(D::LocalGeoData) = embeddim(D.object)
sdata(D::LocalGeoData) = sobj(D) != nothing
nvals(D::LocalGeoData) = nelements(D.object)
snvals(D::LocalGeoData) = nelements(D.data)
spars(D::LocalGeoData) = D.datapars
spars(D::LocalGeoData, i::Int) = spars(D)[i-nvals(D)]
spars(D::LocalGeoData, i::AbstractVector{Int}) = view(spars(D), i .- nvals(D))
nall(D::LocalGeoData) = sdata(D) ? nvals(D)+snvals(D) : nvals(D)
centro(D::LocalGeoData, i::Int) = i<=nvals(D) ? centroid(obj(D),i) : centroid(sobj(D),i-nvals(D))
coords(D::LocalGeoData, i::Int) = coordinates(centro(D,i))
coords(D::SpatialData) = reduce(hcat, [coordinates(centroid(D,x)) for x in 1:nelements(D)])
coords(D::SpatialData, i::AbstractVector{Int}) = reduce(hcat, [coordinates(centroid(D,x)) for x in i])

rotation(D::LocalGeoData) = D.localaniso.rotation
srotation(D::LocalGeoData) = view(rotation(D),spars(D))
magnitude(D::LocalGeoData) = D.localaniso.magnitude
smagnitude(D::LocalGeoData) = view(magnitude(D),:,spars(D))

rotation(D::LocalGeoData, i::Int) = i<=nvals(D) ? rotation(D)[i] : rotation(D)[spars(D,i)]
rotation(D::LocalGeoData, i::AbstractVector{Int}) = i<=nvals(D) ? view(rotation(D), i) : view(rotation(D), spars(D,i))
magnitude(D::LocalGeoData, i::Ints) = i<=nvals(D) ? view(magnitude(D), :, i) : view(magnitude(D), :, spars(D,i))

Base.show(io::IO, ld::LocalGeoData) = print(io, "LocalGeoData")

function Base.show(io::IO, ::MIME"text/plain", ld::LocalGeoData)
	print(io,"LocalGeoData $(ndims(ld))-D")
end

abstract type LocalParMethods end
struct Geometric <: LocalParMethods end
struct Gradients <: LocalParMethods end


abstract type LocalMetric end
struct AnisoDistance  <: LocalMetric end
struct KernelVariogram <: LocalMetric end
struct GraphDistance  <: LocalMetric end
