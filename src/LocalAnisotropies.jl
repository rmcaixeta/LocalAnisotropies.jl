module LocalAnisotropies

using CSV
using Distances
using GeoStatsBase
using ImageFiltering
using KrigingEstimators
using LightGraphs:dijkstra_shortest_paths,connected_components
using LinearAlgebra
using LossFunctions
using Meshes
using MultivariateStats
using NearestNeighbors
using ReferenceFrameRotations
using Setfield
using SimpleWeightedGraphs
using StaticArrays
using StatsBase:Weights,quantile,mean!
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
include("parametrization/conventions.jl")
include("parametrization/conversions.jl")
include("parametrization/geometric.jl")
include("parametrization/gradients.jl")
include("parametrization/interpolation.jl")
include("parametrization/partitions.jl")
include("parametrization/searchers.jl")
include("parametrization/utilities.jl")
include("spacetransforms/deformation.jl")
include("spacetransforms/graph.jl")
include("spacetransforms/metrics.jl")

export
    addgraph,
    convertangles,
    deformspace,
    exportpars,
    localparameters,
    localpars2vtk,
    rescale_magnitude,
    smoothpars,
    Geometric,
    Gradients,
    GraphDistance,
    IDWpars,
    LocalAnisotropy,
    LocalGeoData,
    LocalKriging,
    LocalParameters,
    LocalVariogram,
    TestSet
end
