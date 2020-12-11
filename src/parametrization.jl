
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

struct LocalSpatialData
	object::SpatialData
	localpars::LocalParameters
    graph::SimpleWeightedGraph
end

LocalSpatialData(obj::SpatialData, lpars::LocalParameters) = LocalSpatialData(obj, lpars, nothing)
obj(D::LocalSpatialData) = D.object
rotation(D::LocalSpatialData) = D.localpars.rotation
rotation(D::LocalSpatialData, i::Int) = D.localpars.rotation[i]
rotation(D::LocalSpatialData, i::AbstractVector{Int}) = view(D.localpars.rotation, i)
magnitude(D::LocalSpatialData) = D.localpars.magnitude
magnitude(D::LocalSpatialData, i::Ints) = view(D.localpars.magnitude, :, i)
nvals(D::LocalSpatialData) = nvals(D.object)
ndims(D::LocalSpatialData) = ndims(D.object)

Base.show(io::IO, ld::LocalSpatialData) = print(io, "LocalSpatialData")

function Base.show(io::IO, ::MIME"text/plain", ld::LocalSpatialData)
	println(io,"LocalSpatialData $(ndims(ld))-D")
end


abstract type LocalParametersMethods end
abstract type Geometric <: LocalParametersMethods end
abstract type Gradients <: LocalParametersMethods end


abstract type LocalMetric end
abstract type LocalAnisotropy <: LocalMetric end
abstract type LocalVariogram <: LocalMetric end
abstract type GraphDistance <: LocalMetric end
