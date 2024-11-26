# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsProcesses.jl
# ------------------------------------------------------------------

"""
    LocalSGS(; [paramaters])

The sequential process method introduced by Gomez-Hernandez 1993.
It traverses all locations of the geospatial domain according to a path,
approximates the conditional Gaussian distribution within a neighborhood
using simple Kriging, and assigns a value to the center of the neighborhood
by sampling from this distribution.

## Parameters

* `path`         - Process path (default to `LinearPath()`)
* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `36`)
* `neighborhood` - Search neighborhood (default to `:range`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)
* `init`         - Data initialization method (default to `NearestInit()`)

For each location in the process `path`, a maximum number of
neighbors `maxneighbors` is used to fit the conditional Gaussian
distribution. The neighbors are searched according to a `neighborhood`.

The `neighborhood` can be a `MetricBall`, the symbol `:range` or `nothing`.
The symbol `:range` is converted to `MetricBall(range(γ))` where `γ` is the
variogram of the Gaussian process. If `neighborhood` is `nothing`, the nearest
neighbors are used without additional constraints.

## References

* Gomez-Hernandez & Journel 1993. [Joint Sequential Simulation of
  MultiGaussian Fields](https://link.springer.com/chapter/10.1007/978-94-011-1739-5_8)

### Notes

* This method is very sensitive to the various parameters.
  Care must be taken to make sure that enough neighbors
  are used in the underlying Kriging model.
"""
@kwdef struct LocalSGS{P,N,D,I} <: RandMethod
    localaniso::LocalAnisotropy
    path::P = LinearPath()
    minneighbors::Int = 1
    maxneighbors::Int = 36 # 6x6 grid cells
    neighborhood::N = :range
    distance::D = Euclidean()
    init::I = NearestInit()
end

function GeoStatsProcesses.randprep(
    rng,
    process::GaussianProcess,
    method::LocalSGS,
    setup::RandSetup,
)
    (; path, minneighbors, maxneighbors, neighborhood, distance, init) = method
    method_ = SEQMethod(path, minneighbors, maxneighbors, neighborhood, distance, init)
    GeoStatsProcesses.randprep(rng, process, method_, setup)
end

function GeoStatsProcesses.randsingle(
    rng,
    ::GaussianProcess,
    method::LocalSGS,
    setup::RandSetup,
    prep,
)
    # retrieve parameters
    (; localaniso, path, init) = method
    (; varnames, vartypes) = setup
    (; dom, data, probmodel, marginal, minneighbors, maxneighbors, searcher) = prep

    # initialize buffers for realization and simulation mask
    vars = Dict(zip(varnames, vartypes))
    buff, mask = initbuff(dom, vars, init, data = data)

    # consider point set with centroids for now
    pointset = PointSet([centroid(dom, ind) for ind = 1:nelements(dom)])

    pairs = map(varnames) do var
        # pre-allocate memory for neighbors
        neighbors = Vector{Int}(undef, maxneighbors)

        # retrieve realization and mask for variable
        realization = buff[var]
        simulated = mask[var]

        # simulation loop
        for ind in traverse(dom, path)
            if !simulated[ind]
                center = pointset[ind]
                # search neighbors with simulated data
                nneigh = search!(neighbors, center, searcher, mask = simulated)

                # rebuild probmodel with local par
                localpar = (rotation(localaniso, ind), magnitude(localaniso, ind))
                local_probmodel = mw_estimator(probmodel, localpar)

                if nneigh < minneighbors
                    # draw from marginal
                    realization[ind] = rand(rng, marginal)
                else
                    # neighborhood with data
                    neigh = let
                        ninds = view(neighbors, 1:nneigh)
                        dom = view(pointset, ninds)
                        val = view(realization, ninds)
                        tab = (; var => val)
                        georef(tab, dom)
                    end

                    # fit distribution probmodel
                    fitted = GeoStatsModels.fit(local_probmodel, neigh)

                    # draw from conditional or marginal
                    distribution = if GeoStatsModels.status(fitted)
                        GeoStatsModels.predictprob(fitted, var, center)
                    else
                        marginal
                    end
                    realization[ind] = rand(rng, distribution)
                end

                # mark location as simulated and continue
                simulated[ind] = true
            end
        end

        var => realization
    end

    (; pairs...)
end

# temp solution; to fix in meshes
function Meshes._pboxes(::Type{𝔼{N}}, points) where {N}
    p = first(points)
    ℒ = Meshes.lentype(p)
    cmin = [typemax(ℒ) for i = 1:N]
    cmax = [typemin(ℒ) for i = 1:N]
    for p in points
        c = getfield(coords(p), :coords)
        cmin = [min(c[i], cmin[i]) for i = 1:N]
        cmax = [max(c[i], cmax[i]) for i = 1:N]
    end
    Box(Meshes.withcrs(p, Tuple(cmin)), Meshes.withcrs(p, Tuple(cmax)))
end
