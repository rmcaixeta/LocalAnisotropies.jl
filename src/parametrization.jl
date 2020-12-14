
struct LocalParameters
  rotation::AbstractVector{Quaternion}
  magnitude::AbstractArray{AbstractFloat,2}
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
	println(io,"LocalParameters $(ndims(lp))-D")
end

GeoData = Union{SpatialData,GeoStatsBase.DataView}
struct LocalGeoData
	data::GeoData
	object::GeoData
	localpars::LocalParameters
	refvario::Union{Variogram,Nothing}
	datapars::Union{AbstractVector{Int},Nothing}
    graph::Union{SimpleWeightedGraph,Nothing}
end

function LocalGeoData(obj::GeoData, lpars::LocalParameters)
	LocalGeoData(nothing, obj, lpars, nothing, nothing, nothing)
end

function LocalGeoData(obj::GeoData, lpars::LocalParameters, refvario::Variogram)
	LocalGeoData(nothing, obj, lpars, refvario, nothing, nothing)
end

function LocalGeoData(hd::GeoData, obj::GeoData, lpars::LocalParameters)
	hdids = grid2hd_ids(hd,obj)
	LocalGeoData(hd, obj, lpars, nothing, hdids, nothing)
end

function LocalGeoData(hd::GeoData, obj::GeoData, lpars::LocalParameters, refvario::Variogram)
	hdids = grid2hd_ids(hd,obj)
	LocalGeoData(hd, obj, lpars, refvario, hdids, nothing)
end

obj(D::LocalGeoData) = D.object
sobj(D::LocalGeoData) = D.data
ndims(D::LocalGeoData) = ncoords(D.object)
sdata(D::LocalGeoData) = sobj(D) != nothing
nvals(D::LocalGeoData) = nelms(D.object)
snvals(D::LocalGeoData) = nelms(D.data)
spars(D::LocalGeoData) = D.datapars
spars(D::LocalGeoData, i::Int) = spars(D)[i-nvals(D)]
spars(D::LocalGeoData, i::AbstractVector{Int}) = view(spars(D), i .- nvals(D))
nall(D::LocalGeoData) = sdata(D) ? nvals(D)+snvals(D) : nvals(D)
coords(D::LocalGeoData, i::Int) = i<=nvals(D) ? coordinates(obj(D),i) : coordinates(sobj(D),i-nvals(D))

rotation(D::LocalGeoData) = D.localpars.rotation
srotation(D::LocalGeoData) = view(rotation(D),spars(D))
magnitude(D::LocalGeoData) = D.localpars.magnitude
smagnitude(D::LocalGeoData) = view(magnitude(D),:,spars(D))

rotation(D::LocalGeoData, i::Int) = i<=nvals(D) ? rotation(D)[i] : rotation(D)[spars(D,i)]
rotation(D::LocalGeoData, i::AbstractVector{Int}) = i<=nvals(D) ? view(rotation(D), i) : view(rotation(D), spars(D,i))
magnitude(D::LocalGeoData, i::Ints) = i<=nvals(D) ? view(magnitude(D), :, i) : view(magnitude(D), :, spars(D,i))

Base.show(io::IO, ld::LocalGeoData) = print(io, "LocalGeoData")

function Base.show(io::IO, ::MIME"text/plain", ld::LocalGeoData)
	println(io,"LocalGeoData $(ndims(ld))-D")
end

abstract type LocalParMethods end
struct Geometric <: LocalParMethods end
struct Gradients <: LocalParMethods end


abstract type LocalMetric end
struct LocalAnisotropy <: LocalMetric end
struct LocalVariogram <: LocalMetric end
struct GraphDistance <: LocalMetric end
