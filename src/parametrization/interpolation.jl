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
local neighbors returned from the `searcher`. The interpolation can be inverse
distance weighted by given `power` and `metric` or a simply averaged if
`power=0`. The ellipses/ellipsoids rotations are interpolated using
(weighted) average of quaternions as described by Markley et al (2007). The
interpolated magnitude is the weighted median of the neighbors magnitude.
One `bkgpars` can be assigned where no neighbors are found; must be
informed as (Quaternion, AbstractVector)


## Example

```julia
searcher = KNearestSearch(data, 10)
idw3_into_grid = idwpars(localaniso, searcher, grid, power=3)
```

## Reference

Markley, F.L., et al. (2007). [Averaging quaternions](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872.pdf)
"""
idwpars(lpars, searcher::NeighborSearchMethod, domain; kwargs...) =
    interpolate(lpars, searcher, domain; kwargs...)

"""
    smoothpars(localaniso, searcher; bkgpars=nothing)

Smooth `LocalAnisotropy`, using the local neighbors returned from the `searcher`.
Simple local averages. The smoothed ellipses/ellipsoids rotations are based on
average of quaternions as described by Markley et al. (2007). The interpolated
magnitude is the weighted median of the neighbors magnitude. One `bkgpars`
can be assigned where no neighbors are found; must be informed
as (Quaternion, AbstractVector)

## Example

```julia
searcher = KNearestSearch(data, 10)
averaged_inplace = smoothpars(localaniso, searcher)
```

## Reference

Markley, F.L., et al. (2007). [Averaging quaternions](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070017872.pdf)
"""
smoothpars(lpars, searcher::NeighborSearchMethod; bkgpars = nothing, bkgwgt = 0.0) =
    interpolate(lpars, searcher; power = 0.0, bkgpars, bkgwgt)


function interpolate(
    lpars,
    searcher::NeighborSearchMethod,
    domain = nothing;
    power::Real = 0,
    metric = Euclidean(),
    bkgpars = nothing,
    bkgwgt = 0.0,
)
    D = searcher.domain
    N = embeddim(D)
    smooth = isnothing(domain)
    len = smooth ? nvals(D) : nvals(domain)
    @assert nvals(D) == nvals(lpars) "searcher domain must match number of local anisotropies"

    quat = Array{Quaternion}(undef, len)
    mag = lpars.magnitude
    m = smooth && power == 0 ? mag : Array{Float64}(undef, N, len)

    Threads.@threads for i = 1:len
        ic = smooth ? centro(D, i) : centro(domain, i)
        icoords = ustrip.(to(ic))
        neighids = search(ic, searcher)

        if length(neighids) == 0
            if isnothing(bkgpars)
                throw(
                    ErrorException(
                        "zero neighbors at some location; adjust searcher or set bkgpars",
                    ),
                )
            else
                quat[i] = bkgpars[1]
                m[:, i] .= bkgpars[2]
            end
        elseif smooth
            quat[i] = quatavg(rotation(lpars, neighids))
            #if domain!=nothing
            #  mi      = magnitude(lpars,neighids)
            #  m[:,i] .= mapreduce(x->quantile(view(mi,x,:),0.5), vcat, 1:N)
            #end
        else
            twgt = 1.0 - bkgwgt

            prewgts = if power != 0
                xcoords = coords_(D, neighids)
                1 ./ (eps() .+ Distances.colwise(metric, icoords, xcoords)) .^ power
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

            quat[i] = quatavg(rot, weights)
            wgts = Weights(weights)
            m[:, i] .= mapreduce(x -> quantile(view(mi, x, :), wgts, 0.5), vcat, 1:N)
            # mean instead median for magnitude?
        end
    end

    LocalAnisotropy(quat, m)
end

quat2vector(q) = [q.q0; q.q1; q.q2; q.q3]
tensor(q) = q' .* q
wtensor(q) = q[1] .* (q[2]' .* q[2])

function quatavg(
    qarr::AbstractVector{Quaternion},
    warr::AbstractVector{Float64} = Float64[],
)
    if length(qarr) == 1
        qarr[1]
    else
        Q = mapreduce(quat2vector, hcat, qarr)
        n = size(Q, 2)
        W = length(warr) > 0

        A =
            W ? mapreduce(wtensor, +, zip(warr, eachcol(Q))) :
            mapreduce(tensor, +, eachcol(Q))

        # scale
        Nw = W ? sum(warr) : n
        A ./= Nw

        # compute eigenvalues and -vectors
        T = eigen(Symmetric(SMatrix{4,4}(A)))
        V = view(T.vectors[:, sortperm(T.values, rev = true)], :, 1)
        Quaternion(V...)
    end
end
