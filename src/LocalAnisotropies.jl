module LocalAnisotropies

using Distances
using LinearAlgebra
using ImageFiltering
using MultivariateStats
using NearestNeighbors
using KrigingEstimators
using StaticArrays
using ReferenceFrameRotations
using WriteVTK
using GeoStatsBase
using KrigingEstimators
using Variography

import GeoStatsBase: solve
import KrigingEstimators: FittedKriging, KrigingState, KrigingWeights,
 combine, factorize, nconstraints, set_constraints_rhs!, set_constraints_lhs!

include("estimators.jl")
include("parametrization.jl")
include("partitions.jl")
include("variograms.jl")

export
    LocalKriging,
    gradients

end
