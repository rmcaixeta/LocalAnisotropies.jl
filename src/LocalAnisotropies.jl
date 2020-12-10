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
include("estimators/cross_validation.jl")
include("estimators/idw.jl")
include("estimators/kriging.jl")
include("estimators/variograms.jl")
include("parametrization/calibration.jl")
include("parametrization/geometric.jl")
include("parametrization/gradients.jl")
include("parametrization/interpolation.jl")
include("parametrization/partitions.jl")
include("parametrization/searchers.jl")
include("parametrization/utilities.jl")
#include("spatialtransforms/lmds.jl")
#include("spatialtransforms/metrics.jl")

export
    localpars,
    localpars2vtk,
    rescale_magnitude,
    smooth,
    Geometric,
    Gradients,
    LocalKriging,
    TestSet
end
