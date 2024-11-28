# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module LocalAnisotropies

using Distances
using GeoStatsBase
using GeoStatsFunctions
using GeoStatsModels
using GeoStatsProcesses
using GeoTables
using Graphs: dijkstra_shortest_paths, connected_components
using ImageFiltering
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
using WriteVTK

import Distributions: Normal
import GeoStatsModels:
    FittedIDW,
    FittedKriging,
    GeoStatsModel,
    IDW,
    IDWState,
    KrigingModel,
    KrigingState,
    KrigingWeights,
    OrdinaryKriging,
    SimpleKriging,
    lhs,
    nconstraints,
    predict,
    predictmean,
    predictprob,
    predictvar,
    set_constraints_rhs!
import GeoStatsProcesses: RandMethod, RandSetup
import GeoStatsTransforms: ColumnSelector, TableTransform, selector, apply

include("parametrization.jl")
include("estimators/idw.jl")
include("estimators/kriging.jl")
include("estimators.jl")
include("estimators/variograms.jl")
include("estimators/sgs.jl")
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
    KernelVariogram,
    LocalAnisotropy,
    LocalGeoData,
    LocalIDW,
    LocalInterpolate,
    LocalKriging,
    LocalSGS,
    RotationRule
end
