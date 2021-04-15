# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    LocalParameters(rotation, magnitude)

Create a `LocalParameters` object. The `rotation` is a vector of
`ReferenceFrameRotations.Quaternion` and the magnitude is `AbstractArray{Float}`
of size (number of dimensions, number of elements)
"""
struct LocalParameters
  rotation::AbstractVector{Quaternion}
  magnitude::AbstractArray
end

Ints = Union{Int,AbstractVector{Int}}
rotation(L::LocalParameters) = L.rotation
rotation(L::LocalParameters, i::Int) = L.rotation[i]
rotation(L::LocalParameters, i::AbstractVector{Int}) = view(L.rotation, i)
magnitude(L::LocalParameters) = L.magnitude
magnitude(L::LocalParameters, i::Ints) = view(L.magnitude, :, i)
slice(L::LocalParameters, i::Ints) = LocalParameters(rotation(L,i),magnitude(L,i))
nvals(L::LocalParameters) = length(L.rotation)
ndims(L::LocalParameters) = size(L.magnitude,1)

Base.show(io::IO, lp::LocalParameters) = print(io, "LocalParameters")

function Base.show(io::IO, ::MIME"text/plain", lp::LocalParameters)
	print(io,"LocalParameters $(ndims(lp))-D")
end

"""
    SpatialData

Union of multiple types of spatial data from GeoStats.jl: `GeoData`, `Domain`,
`DataView` and `DomainView`
"""
SpatialData = Union{GeoData,Domain,DataView,DomainView}

struct LocalGeoData
	data::SpatialData
	object::SpatialData
	localpars::LocalParameters
	refvario::Union{Variogram,Nothing}
	datapars::Union{AbstractVector{Int},Nothing}
    graph::Union{SimpleWeightedGraph,Nothing}
end

"""
    LocalGeoData(object, localpars)

Association of a domain object and its `LocalParameters`.
"""
function LocalGeoData(obj::SpatialData, lpars::LocalParameters)
	LocalGeoData(nothing, obj, lpars, nothing, nothing, nothing)
end

"""
    LocalGeoData(object, localpars, refvario)

Association of a domain object, its `LocalParameters` and its reference variogram.
"""
function LocalGeoData(obj::SpatialData, lpars::LocalParameters, refvario::Variogram)
	LocalGeoData(nothing, obj, lpars, refvario, nothing, nothing)
end

"""
    LocalGeoData(samples, object, localpars)

Association of samples, a domain object and their `LocalParameters`. The samples
`LocalParameters` are passed from the domain via nearest neighbors.
"""
function LocalGeoData(hd::SpatialData, obj::SpatialData, lpars::LocalParameters)
	hdids = grid2hd_ids(hd,obj)
	LocalGeoData(hd, obj, lpars, nothing, hdids, nothing)
end

"""
    LocalGeoData(samples, object, localpars, refvario)

Association of samples, a domain object, their `LocalParameters` and the
reference variogram. The samples `LocalParameters` are passed from the domain
via nearest neighbors.
"""
function LocalGeoData(hd::SpatialData, obj::SpatialData, lpars::LocalParameters, refvario::Variogram)
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

rotation(D::LocalGeoData) = D.localpars.rotation
srotation(D::LocalGeoData) = view(rotation(D),spars(D))
magnitude(D::LocalGeoData) = D.localpars.magnitude
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
struct LocalAnisotropy <: LocalMetric end
struct LocalVariogram  <: LocalMetric end
struct GraphDistance   <: LocalMetric end
