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
    targetD = smooth ? D : domain
    #@assert nvals(D) == nvals(lpars) "searcher domain must match number of local anisotropies"

    quat = Array{Quaternion}(undef, len)
    mag = lpars.magnitude
    m = Array{Float64}(undef, N, len)

    Threads.@threads for i = 1:len
        ic = centro(targetD, i)
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
        else
            twgt = 1.0 - bkgwgt

            # function to return weights according to method
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

            #quat[i], mx = ellipsavg(mi, rot, weights)
            #m[:,i] .= mx
            #m[2,i] = m[1,i] # if simplify ....

            quat[i], f = quatavg(rot, weights)
            wgts = Weights(weights)
            mx = mapreduce(x -> quantile(view(mi, x, :), wgts, 0.5), vcat, 1:N)
            mx = f .* mx .+ (1 - f) .* [1.0 for i in 1:N]
            m[:, i] .= mx
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
        qarr[1], 1.0
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
        s = sortperm(T.values, rev = true)
        f = T.values[s][1] / sum(T.values)
        V = view(T.vectors[:, s], :, 1)
        Quaternion(V...), f
    end
end

function ellipsavg(
    mag::AbstractArray,
    qarr::AbstractVector{Quaternion},
    wgt::AbstractVector{Float64} = Float64[],
)
    nv = length(qarr)
    if nv == 1
        qarr[1], mag
    else
        Ms = [LocalAnisotropies.qmat(qarr[i], mag[:, i]) for i = 1:nv]
        wgt = length(wgt) == 0 ? [1 / nv for i in Ms] : wgt
        log_Ms = [log(Matrix(M)) for M in Ms]
        log_mean = sum(wgt[i] * log_Ms[i] for i = 1:length(Ms))
        M_avg = exp(log_mean)

        N = size(M_avg, 1)
        T = eigen(Symmetric(SMatrix{N,N}(M_avg)))
        s = sortperm(T.values, rev = true)
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
