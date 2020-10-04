module LocalAnisotropies

using Distances
using LinearAlgebra
using MultivariateStats
using NearestNeighbors
using GeoStats
using StaticArrays
using ReferenceFrameRotations
using WriteVTK

include("estimators.jl")
include("parametrization.jl")
include("partitions.jl")
include("variograms.jl")

#export
#    function1,
#    function1

end
