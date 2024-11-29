# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsProcesses.jl
# ------------------------------------------------------------------

"""
    LocalSGS(; [paramaters])
    LocalSGS(localaniso; [paramaters])
    LocalSGS(method, localaniso; [paramaters])

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
    method::Symbol = :MovingWindows
    localaniso::LocalAnisotropy
    path::P = LinearPath()
    minneighbors::Int = 1
    maxneighbors::Int = 36 # 6x6 grid cells
    neighborhood::N = :range
    distance::D = Euclidean()
    init::I = NearestInit()
end

LocalSGS(localaniso::LocalAnisotropy; kwargs...) = LocalSGS(; localaniso, kwargs...)
LocalSGS(method::Symbol, localaniso::LocalAnisotropy; kwargs...) =
    LocalSGS(; method, localaniso, kwargs...)

function GeoStatsProcesses.randprep(
    rng,
    process::GaussianProcess,
    meth::LocalSGS,
    setup::RandSetup,
)
    (; path, minneighbors, maxneighbors, neighborhood, distance, init) = meth
    method_ = SEQMethod(path, minneighbors, maxneighbors, neighborhood, distance, init)
    GeoStatsProcesses.randprep(rng, process, method_, setup)
end

function GeoStatsProcesses.randsingle(
    rng,
    ::GaussianProcess,
    meth::LocalSGS,
    setup::RandSetup,
    prep,
)
    # retrieve parameters
    (; method, localaniso, path, init) = meth
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

                    # rebuild probmodel as a local probmodel
                    hdlocalaniso = neighs_localaniso(localaniso, dom, neigh; method)
                    local_model = local_probmodel(probmodel, localaniso, hdlocalaniso)

                    # fit distribution probmodel
                    fitted = local_fit(local_model, neigh, i = ind, m = 1:nrow(neigh))

                    # draw from conditional or marginal
                    distribution = if local_status(fitted)
                        predictprob(fitted, var, center)
                    else
                        @warn "Kriging error at point $ind; drawing random value"
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

function local_probmodel(probmodel, localaniso, hdlocalaniso)
    isnothing(hdlocalaniso) ? MW_SKModel(localaniso, probmodel.γ, probmodel.μ) :
        KC_SKModel(localaniso, probmodel.γ, probmodel.μ, hdlocalaniso)
end

function Meshes._pboxes(::Type{𝔼{N}}, points) where {N}
    p = first(points)
    ℒ = lentype(p)
    cmin = fill(typemax(ℒ), N)
    cmax = fill(typemin(ℒ), N)
  
    for p in points
      c = getfield(coords(p), :coords)
      for i in 1:N
        cmin[i] = min(c[i], cmin[i])
        cmax[i] = max(c[i], cmax[i])
      end
    end
    Box(withcrs(p, Tuple(cmin)), withcrs(p, Tuple(cmax)))
  end