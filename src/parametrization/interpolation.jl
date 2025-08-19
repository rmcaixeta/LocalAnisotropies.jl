# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    nnpars(lpars, data, target)

Pass local anisotropies `lpars` from `data` to `target` object, via nearest
neighbors. Useful to pass parameters from domain to samples.
"""
nnpars(lpars, data, target) = slice(lpars, grid2hd_ids(target, data))

"""
    idwpars(localaniso, searcher, domain; power=2, metric=Euclidean())

Interpolate `LocalAnisotropy` data into a `domain`, using the
local neighbors returned from the `searcher`. Weights calculated via inverse
distance weighted by given `power` and `metric`. Custom kwargs explained in
`interpolate`


## Example

```julia
searcher = KNearestSearch(data, 10)
idw3_into_grid = idwpars(localaniso, searcher, grid, power=3)
```
"""
idwpars(lpars, searcher::NeighborSearchMethod, domain; kwargs...) =
  interpolate(lpars, searcher, domain; wgtmethod=:idw, kwargs...)

"""
    smoothpars(localaniso, searcher; b=0.0, kwargs)
    smoothpars(localaniso, searcher, domain; b=0.0, kwargs)

Smooth `LocalAnisotropy`, using the local neighbors returned from the `searcher`.
Simple local averages if `b` = 0 otherwise gaussian kernel smoothing with given
bandwidth. Custom kwargs explained in `interpolate`

## Example

```julia
searcher = KNearestSearch(data, 10)
averaged_inplace = smoothpars(localaniso, searcher)
averaged_grid = smoothpars(localaniso, searcher, grid, b=0.6)
```
"""
smoothpars(lpars, searcher::NeighborSearchMethod, domain=nothing; b=0.0, kwargs...) =
  interpolate(lpars, searcher, domain; wgtmethod=:kernel, power=b, kwargs...)

"""
    interpolate(localaniso, searcher, domain=nothing;
        wgtmethod = :idw, power=0, metric=Euclidean(),
        bkgpars=nothing, bkgwgt=0.0, fillna=:bkg, method=:qavg)

General function to return the (weighted) average of `localaniso` using `searcher`.

## Keyword Parameters

* `wgtmethod` - :idw for inverse distance weighting or :kernel for gaussian kernel smoothing
* `power`     - power for :idw or the bandwidth for :kernel
* `metric`    - metric for distances calculation
* `bkgpars`   - default rotation and magnitude informed as (Quaternion, AbstractVector)
* `bkgwgt`    - weight between [0,1] assigned to default background parameters
* `fillna`    - :bkg will assign the `bkgpars` if searcher returns nothing and :nn
  will run a nearest neighbor with what was estimated successfully - also used
  if `bkgpars` is not informed
* `method`    - how the average of rotation and magnitude is done; options below:

  - :qavg, :qavg_mn, :qavg_md, :qavg_full, :qavg_mn_full, :qavg_md_full

  The ellipses/ellipsoids rotations are interpolated using (weighted) average of
  quaternions as described by Markley et al (2007). Magnitude is kept equal if
  :qavg and the data is not interpolated to a domain. The sufix "_mn" refers to
  weighted mean of the neighbors magnitude and "_md" to the median. The sufix "_full"
  will add an isotropic component in the magnitude interpolation, the more the wights to it
  the more the rotations lack similiarity.

  - :ellipavg, :ellipavg_full

  The ellipses/ellipsoids are converted to rotation matrix scaled by the magnitude
  and a Log-Euclidean (Weighted) Mean is applied followed by extracting eigen-
  values and vectors. Eigenvectors will represent the new rotation and eigenvalues the
  magnitude. The only difference between :ellipavg and :ellipavg_full is that :ellipavg
  will return simplified magnitude for 3-D cases, making `r1`=`r2`

## Reference

Markley, F.L., et al. (2007). [Averaging quaternions](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872.pdf)
"""

function interpolate(
  lpars,
  searcher::NeighborSearchMethod,
  domain=nothing;
  wgtmethod::Symbol=:idw,
  power::Real=0,
  method::Symbol=:qavg,
  metric=Euclidean(),
  bkgpars=nothing,
  bkgwgt=0.0,
  fillna=:bkg
)
  D = searcher.domain
  N = embeddim(D)
  smooth = isnothing(domain)
  len = smooth ? nvals(D) : nvals(domain)
  targetD = smooth ? D : domain
  #@assert nvals(D) == nvals(lpars) "searcher domain must match number of local anisotropies"

  quat = Vector{Quaternion}(undef, len)
  m = Array{Float64}(undef, N, len)
  missids = zeros(Bool, len)

  @tasks for i in 1:len
    ic = centro(targetD, i)
    icoords = ustrip.(to(ic))
    neighids = search(ic, searcher)

    if length(neighids) == 0
      missids[i] = 1
    else
      twgt = 1.0 - bkgwgt

      prewgts = if power != 0
        xcoords = coords_(D, neighids)
        distances = Distances.colwise(metric, icoords, xcoords)
        get_weights(wgtmethod, distances, power)
      else
        [1.0 for i in neighids]
      end
      weights = twgt .* prewgts ./ sum(prewgts)

      rot = rotation(lpars, neighids)
      mi = magnitude(lpars, neighids)
      if bkgwgt > 0
        rot = vcat(bkgpars[1], rot)
        mi = hcat(bkgpars[2], mi)
        weights = vcat(bkgwgt, weights)
      end

      if method in (:ellipavg, :ellipavg_full)
        quat[i], mx = ellipsavg(mi, rot, weights)
        m[:, i] .= mx
        if method == :ellipavg && N == 3
          newm = mean(m[1:2, i])
          m[1, i] = 1
          m[2, i] = 1
          m[3, i] = m[3, i] / newm
        end
      else
        quat[i], f = occursin("_md", "$method") ? quatmed(rot, weights) : quatavg(rot, weights)
        if smooth && method == :qavg
          m[:, i] .= magnitude(lpars, i)
        else
          wgts = Weights(weights)
          mx =
            occursin("_mn", "$method") ? (mapreduce(x -> mean(view(mi, x, :), wgts), vcat, 1:N)) :
            (mapreduce(x -> quantile(view(mi, x, :), wgts, 0.5), vcat, 1:N))
          mx = occursin("_full", "$method") ? (f .* mx .+ (1 - f) .* [1.0 for i in 1:N]) : mx
          m[:, i] .= mx
        end
      end
    end
  end

  # deal with missing
  ismiss = findall(missids)
  if length(ismiss) > 0
    if fillna == :nn || isnothing(bkgpars)
      notmiss = findall(.!missids)
      Dfrom = view(targetD, notmiss)
      Dto = view(targetD, ismiss)
      outids = notmiss[grid2hd_ids(Dto, Dfrom)]
      quat[ismiss] .= quat[outids]
      m[:, ismiss] .= m[:, outids]
    else
      for i in ismiss
        quat[i] = bkgpars[1]
        m[:, i] .= bkgpars[2]
      end
    end
  end

  LocalAnisotropy(quat, m)
end

quat2vector(q) = [q.q0; q.q1; q.q2; q.q3]
tensor(q) = q' .* q
wtensor(q) = q[1] .* (q[2]' .* q[2])

function quatavg(qarr::AbstractVector{Quaternion}, warr::AbstractVector{Float64}=Float64[])
  if length(qarr) == 1
    qarr[1], 1.0
  else
    Q = mapreduce(quat2vector, hcat, qarr)
    n = size(Q, 2)
    W = length(warr) > 0

    A = W ? mapreduce(wtensor, +, zip(warr, eachcol(Q))) : mapreduce(tensor, +, eachcol(Q))

    # scale
    Nw = W ? sum(warr) : n
    A ./= Nw

    # compute eigenvalues and -vectors
    T = eigen(Symmetric(SMatrix{4,4}(A)))
    s = sortperm(T.values, rev=true)
    f = T.values[s][1] / sum(T.values)
    V = view(T.vectors[:, s], :, 1)
    Quaternion(V...), f
  end
end

function quatmed(qarr::AbstractVector{Quaternion}, warr::AbstractVector{Float64}=Float64[])
  if length(qarr) == 1
    qarr[1], 1.0
  else
    n = mapreduce(vcat, qarr) do q
      to_vector(q, 3)
    end
    i_idx = argmin(abs.(n[:,1] .- median(n[:,1])))
    j_idx = argmin(abs.(n[:,2] .- median(n[:,2])))
    k_idx = argmin(abs.(n[:,3] .- median(n[:,3])))

    i_idx = findall(n[:,1] .== n[i_idx,1])
    j_idx = findall(n[:,2] .== n[j_idx,2])
    k_idx = findall(n[:,3] .== n[k_idx,3])

    idxs = unique(vcat(i_idx,j_idx,k_idx))
    wgts = isempty(warr) ? warr : view(warr,idxs)
    quatavg(view(qarr,idxs), wgts)
  end
end

function ellipsavg(mag::AbstractArray, qarr::AbstractVector{Quaternion}, wgt::AbstractVector{Float64}=Float64[])
  nv = length(qarr)
  if nv == 1
    qarr[1], mag
  else
    Ms = [qmat(qarr[i], mag[:, i]) for i in 1:nv]
    wgt = length(wgt) == 0 ? [1 / nv for i in Ms] : wgt
    log_Ms = [log(Matrix(M)) for M in Ms]
    log_mean = sum(wgt[i] * log_Ms[i] for i in 1:length(Ms))
    M_avg = exp(log_mean)

    N = size(M_avg, 1)
    T = eigen(Symmetric(SMatrix{N,N}(M_avg)))
    s = sortperm(T.values, rev=true)
    V = T.vectors[:, sortperm(T.values)]'
    λ = T.values[s] .^ 0.5
    λ = λ / λ[1]

    if N == 3
      eigv = V
      det(V) < 0 && (eigv = Diagonal(SVector{3}([-1, 1, 1])) * eigv)
    else
      eigv = zeros(Float64, 3, 3)
      eigv[1:2, 1:2] = V
      eigv[1, 1] ≉ eigv[2, 2] && (eigv[1, :] .*= -1)
      eigv[3, 3] = 1.0
      eigv = SMatrix{3,3}(eigv)
    end
    q = dcm_to_quat(DCM(eigv))
    q, λ
  end
end

function get_weights(method, distances, par)
  method == :kernel ? exp.(-distances .^ 2 ./ (2 * par^2)) :
  method == :idw ? 1 ./ (eps() .+ distances) .^ par : 1 ./ (eps() .+ distances) .^ par
end
