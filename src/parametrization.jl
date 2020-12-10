"""
 Extract local anisotropy from some reference input
 - Gradients 2D and 3D
 - PCA for scattered points

 Store parameters as quaternions + ratios

- Need to resolve:
- best design to scale ratios easily?
    fix reference vario and scale between ]0,1]?
    fix vario range to one and scale between [min,max] ranges?
- how to assign local pars differently to each structures? mayber better to only
allow if the ratios should rescale all or given list of structures


"""

struct LocalParameters
  rotation::AbstractVector{Quaternion}
  magnitude::AbstractArray{AbstractFloat,2}
end

Base.show(io::IO, dh::LocalParameters) = print(io, "LocalParameters")

function Base.show(io::IO, ::MIME"text/plain", dh::LocalParameters)
	N, len = size(dh.magnitude)
	println(io,"LocalParameters $(N)-D")
end

slicelp(lp::LocalParameters,ids) = LocalParameters(lp.rotation[ids],lp.magnitude[:,ids])
