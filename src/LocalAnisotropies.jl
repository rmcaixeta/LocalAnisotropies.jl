# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module LocalAnisotropies

using CSV
using Distances
using GeoStatsBase
using GeoTables
using ImageFiltering
using GeoStatsModels
using Graphs:dijkstra_shortest_paths,connected_components
using LinearAlgebra
using Meshes
using MultivariateStats
using NearestNeighbors
using RecipesBase
using ReferenceFrameRotations
using Setfield
using SimpleWeightedGraphs
using StaticArrays
using StatsBase:Weights,quantile,mean!
using Tables
using Variography
using WriteVTK

import GeoStatsBase: solve
import GeoStatsModels: FittedKriging, KrigingState, KrigingWeights,
 predictmean, predictvar, factorize, nconstraints, set_constraints_rhs!, set_constraints_lhs!

include("estimators.jl")
include("parametrization.jl")
include("estimators/idw.jl")
include("estimators/kriging.jl")
include("estimators/variograms.jl")
include("parametrization/conventions.jl")
include("parametrization/conversions.jl")
include("parametrization/geometric.jl")
include("parametrization/gradients.jl")
include("parametrization/interpolation.jl")
include("parametrization/recipes.jl")
include("parametrization/searchers.jl")
include("parametrization/utilities.jl")
include("parametrization/variograms.jl")
include("spacetransforms/deformation.jl")
include("spacetransforms/graph.jl")
include("spacetransforms/metrics.jl")

export
    convertangles,
    deformspace,
    exportpars,
    graph,
    idwpars,
    localanisotropies,
    localvariography,
    magnitude,
    nnpars,
    reference_magnitude, reference_magnitude!,
    rescale_magnitude, rescale_magnitude!,
    smoothpars,
    to_3d,
    to_vtk,
    AnisoDistance,
    Geometrical,
    Gradients,
    GraphDistance,
    LocalAnisotropy,
    LocalGeoData,
    LocalKriging,
    KernelVariogram,
    RotationRule
end
