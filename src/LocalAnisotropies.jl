module LocalAnisotropies

using Distances
using GeoStatsBase
using ImageFiltering
using KrigingEstimators
using KrigingEstimators
using LinearAlgebra
using LossFunctions
using MultivariateStats
using NearestNeighbors
using ReferenceFrameRotations
using Setfield
using StaticArrays
using Variography
using WriteVTK

import GeoStatsBase: solve
import KrigingEstimators: FittedKriging, KrigingState, KrigingWeights,
 combine, factorize, nconstraints, set_constraints_rhs!, set_constraints_lhs!

include("estimators.jl")
include("parametrization.jl")
include("partitions.jl")
include("quaternions.jl")
include("variograms.jl")
include("tmp/cross_validation.jl")

export
    gradients,
    localpars,
    localpars2vtk,
    pcavector,
    rescale_magnitude,
    smooth,
    LocalKriging,
    TestSet
end
