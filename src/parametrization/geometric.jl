# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    localanisotropies(Geometric, searcher; simplify=true)

Extract `LocalAnisotropy` from a searcher object. Based on the spatial data and
the parameters from the searcher object, a group of neighbor points is collected
at each place. PCA is applied to their spatial coordinates, returning eigenvectors
that represent a local ellipse/ellipsoid aligned with the directions with more
clustered points. The magnitude is the ratio of the PCA eigenvalues. It's useful
to extract local anisotropies from 3-D surface points. If `simplify = true` and
data is 3-D, only the minimum direction is extracted from PCA; the other axes
are calculated (main direction will be always along max dip direction and the
intermediate axis will be orthogonal to the others). This is default for 3-D.

## Example

```julia
searcher = KNearestSearch(pts3dsurface, 10)
prelpars = localanisotropies(Geometric, searcher) # extract local pars at surface
lpars = idwpars(prelpars, searcher, grid) # interpolate local pars to a grid
```
"""
function localanisotropies(
    ::Type{Geometric},
    searcher::NeighborSearchMethod;
    simplify::Bool = true,
)
    D = searcher.domain
    X = coords_(D)
    N, len = size(X)

    quat = Array{Quaternion}(undef, len)
    m = Array{Vector}(undef, len)

    @tasks for i = 1:len
        neighids = search(centro(D, i), searcher)
        λ, v = pca(view(X, :, neighids), simplify)

        det(v) < 0 && (v = Diagonal(SVector{3}([-1, 1, 1])) * v)

        q = dcm_to_quat(DCM(v))
        quat[i] = q
        m[i] = λ / λ[1]
    end

    # deal with -1 eigvals
    m = reduce(hcat, vec(m))
    posm = m .> 0
    if sum(posm) < N * len
        for d = 2:N
            posd = view(posm, d, :)
            sum(posd) == len && continue
            minx = minimum(view(m, d, findall(posd)))
            m[d, findall(.!posd)] .= minx
        end
    end

    LocalAnisotropy(quat, m)
end

function localanisotropies(::Type{Geometric}, trs::GeometrySet, magnitude::AbstractVector)
    q = mapreduce(vcat, trs) do tr
        v = vertices(tr)
        c = centroid(tr)
        n = [ustrip(i) for i in normal(tr)]
        p = minpt(tr)
        v1 = [ustrip(i) for i in (p - c)]
        v1 ./= sum(v1)
        v2 = cross(v1, n)

        m = SMatrix{3,3}(v1..., v2..., n...)
        det(m) < 0 && (m = Diagonal(SVector{3}([-1, 1, 1])) * m)

        dcm_to_quat(DCM(m))
    end
    mag = repeat(magnitude, 1, length(q))
    LocalAnisotropy(q, mag)
end

function pca(X, simplify)
    N = size(X, 1)
    M = MultivariateStats.fit(PCA, ustrip.(X), maxoutdim = N, pratio = 1)
    λ = principalvars(M) ./ principalvars(M)[1]
    nv = length(λ)
    v = N == 3 ? projection(M) : vcat(projection(M), [0 0 1][1:nv])

    if nv == 1
        vx = [-v[2, 1]; v[1, 1]; v[3, 1]]
        v = hcat(v, vx, cross(v[:, 1], vx))
        append!(λ, [λ[1], λ[1]])
    elseif nv == 2 && N == 3
        v = hcat(v, cross(v[:, 1], v[:, 2]))
        push!(λ, minimum(λ))
    end

    if simplify && N == 3
        # use maxdip vector and cross(maxdip,normal) as the main vectors
        abs(v[3, 3]) > 1 && (v[3, 3] = round(v[3, 3], digits = 1))
        az, dp = [atan(v[1, 3], v[2, 3]), -asin(v[3, 3])]
        if dp < 0
            dp += pi / 2
        else
            dp -= pi / 2
        end
        v[:, 1] .= [sin(az) * cos(dp), cos(az) * cos(dp), -sin(dp)]
        v[:, 2] .= cross(v[:, 1], v[:, 3])
        λ[1:2] .= [1.0, 1.0]
    end

    λ[1:N], SMatrix{3,3}(v')
end


function minpt(tr)
    minz = to(boundingbox(tr).min)[3]
    vc = vertices(tr)
    m = [to(p)[3] == minz for p in vc]
    out = vc[m]
    out isa Point ? out : out[1]
end
