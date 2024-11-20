# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    localvariography(sdata, lpars, var₁, var₂=var₁;
    axis=:X, tol=1e-6, maxratio1=Inf, maxratio2=Inf, [optional parameters])

Computes an unconventional directional (cross-)variogram for the variables `var₁`
and `var₂` of spatial data `sdata`. This uses local anisotropies `lpars` to create
pseudo-partitions along local `axis` direction (`axis` ∈ [:X, :Y, :Z]) and the
EmpiricalVariogram is calculated with these data. The directional tolerance
bandwidth `tol` can be passed. The `maxratio1` and `maxratio2` partially ignore
data if magnitude ratios are above these values (does not ignore anything by
default). Other optional parameters are defined at `GeoStats.EmpiricalVariogram`,
such as `nlags`, `maxlag`, `distance` and `algo`. If axis informed is a combination
of two axes (e.g. :XY), this will consider a planar variogram along that plane.

Similar (but not equal) to https://github.com/rmcaixeta/Local_variography
"""
function localvariography(
    obj::SpatialData,
    lpars::LocalAnisotropy,
    var₁,
    var₂ = var₁;
    axis = :X,
    tol = 1e-6,
    maxratio1 = Inf,
    maxratio2 = Inf,
    kwargs...,
)
    p = pseudolocalpartition(obj, lpars, axis, tol, maxratio1, maxratio2)
    EmpiricalVariogram(p, var₁, var₂; kwargs...)
end

function pseudolocalpartition(obj, lpars, axis, tol, maxratio1, maxratio2)
    subs = []
    dims = ndims(lpars)
    n = nvals(obj)
    planar = length("$axis") == 2
    for i = 1:n
        maxratio1 != Inf && ratio1(lpars, i) > maxratio1 && continue
        maxratio2 != Inf && ratio2(lpars, i) > maxratio2 && continue
        x = to(centro(obj, i))
        axfun = planar ? get_normal_axis : iaxis
        v = rotmat(lpars, i)[axfun(axis), 1:dims]
        ptfun = planar ? PlanePartition : DirectionPartition
        p = ptfun(Tuple(v), tol = tol)
        s = Int[i]
        for j = 1:n
            i == j && continue
            y = to(centro(obj, j))
            p(x, y) && push!(s, j)
        end
        push!(subs, s)
    end

    Partition(obj, subs)
end

function get_normal_axis(plan)
    all_ax = Set(['X', 'Y', 'Z'])
    plan_ax = Set("$plan")
    normal_ax = setdiff(all_ax, plan_ax)
    iaxis(Symbol(first(normal_ax)))
end
