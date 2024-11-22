# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module LocalAnisotropies

using Distances
using GeoStatsBase
using GeoTables
using ImageFiltering
using GeoStatsModels
using Graphs: dijkstra_shortest_paths, connected_components
using LinearAlgebra
using Meshes
using MultivariateStats
using NearestNeighbors
using OhMyThreads: tmap, @tasks
using RecipesBase
using ReferenceFrameRotations
using Setfield
using SimpleWeightedGraphs
using StaticArrays
using StatsBase: Weights, quantile, mean!
using Tables
using Unitful: ustrip
using GeoStatsFunctions
using WriteVTK

import GeoStatsModels:
    FittedKriging,
    IDW,
    IDWState,
    FittedIDW,
    KrigingState,
    KrigingWeights,
    KrigingModel,
    GeoStatsModel,
    lhs,
    predict,
    nconstraints,
    predictmean,
    predictvar,
    set_constraints_rhs!
import GeoStatsTransforms: ColumnSelector, TableTransform, selector, apply

include("parametrization.jl")
include("estimators/idw.jl")
include("estimators/kriging.jl")
include("estimators.jl")
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

export adjust_rake,
    adjust_rake!,
    convertangles,
    deformspace,
    graph,
    idwpars,
    localanisotropies,
    localvariography,
    magnitude,
    nnpars,
    reference_magnitude,
    reference_magnitude!,
    rescale_magnitude,
    rescale_magnitude!,
    smoothpars,
    to_3d,
    to_table,
    to_vtk,
    AnisoDistance,
    Geometric,
    Gradients,
    GraphDistance,
    LocalInterpolate,
    LocalAnisotropy,
    LocalGeoData,
    LocalKriging,
    LocalIDW,
    KernelVariogram,
    RotationRule
end
