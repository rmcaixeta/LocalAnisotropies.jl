# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function colwise(D::LocalGeoData, ::Type{AnisoDistance}, i::Int, xj::Vector{Int})
  Qi, ix = qmat(D, i), coords_(D, i)
  d = Vector{Float64}()
  for j in xj
    Qj = qmat(D, j)
    Qij = (Qi + Qj) / 2
    jx = coords_(D, j)
    push!(d, Distances.evaluate(Mahalanobis(Symmetric(Qij)), ix, jx))
  end
  d
end

function colwise(D::LocalGeoData, ::Type{KernelVariogram}, i::Int, xj::Vector{Int})
  Qi, ix = qmat(D, i), centro(D, i)
  d = Vector{Float64}()
  for j in xj
    Qj, jx = qmat(D, j), centro(D, j)
    dγ = abs(sill(D.refvario) - kccov(D.refvario, ix, jx, Qi, Qj))
    push!(d, dγ)
  end
  d
end

function colwise(D::LocalGeoData, ::Type{GraphDistance}, i::Int, j::Ints)
  d = dijkstra_shortest_paths(D.graph, i).dists
  view(d, j)
end

function evaluate(D::LocalGeoData, ::Type{AnisoDistance}, i::Int, j::Int)
  Qi, Qj = [qmat(D, x) for x in (i, j)]
  ix, jx = [coords_(D, x) for x in (i, j)]
  Qij = (Qi + Qj) / 2
  Distances.evaluate(Mahalanobis(Symmetric(Qij)), ix, jx)
end

function evaluate(D::LocalGeoData, ::Type{KernelVariogram}, i::Int, j::Int)
  Qi, Qj = [qmat(D, x) for x in (i, j)]
  ix, jx = [centro(D, x) for x in (i, j)]
  abs(sill(D.refvario) - kccov(D.refvario, ix, jx, Qi, Qj))
end

evaluate(D::LocalGeoData, ::Type{GraphDistance}, i::Int, j::Int) = colwise(GraphDistance, D, i, j)[1]
