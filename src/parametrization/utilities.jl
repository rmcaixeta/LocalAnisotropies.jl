# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    rescale_magnitude(localaniso; r1, r2, r3, clip=[0.05,0.95])
    rescale_magnitude!(localaniso; r1, r2, r3, clip=[0.05,0.95])

Rescale magnitude values of `localaniso` to desired limits using min-max scaling.
`r1`, `r2` and `r3` refer to the ranges along each main axis. `clip` defines
the quantiles to which min and max limits are defined. The default is set to
[0.05,0.95] to avoid outliers influence.

## Example

In the example below, a 3-D local anisotropy object is rescaled. `r1` has all the
values set to 1.0. `r2` and `r3` are rescaled between (0.2, 1.0) and (0.1, 0.5)
respectively. That means that if we associate it to a SphericalVariogram of range=10,
range1 will be equal to 10 at every place. range2 will vary between 2 and 10.
And range3 will vary between 1 and 5.

```julia
example = rescale_magnitude(localaniso3d, r1=1, r2=(0.2,1.0), r3=(0.1,0.5))
```
"""
function rescale_magnitude!(
    lp::LocalAnisotropy;
    r1 = nothing,
    r2 = nothing,
    r3 = nothing,
    clip = [0.05, 0.95],
)
    N = ndims(lp)
    m = lp.magnitude
    rvals = (r1, r2, r3)

    for i = 1:N
        rvals[i] == nothing && continue

        mi = m[i, :]
        qx = quantile(mi, clip)

        if rvals[i] isa Number
            m[i, :] .= rvals[i]
        elseif qx[1] == qx[2]
            throw(ArgumentError("clipped values are unique; inform a single value to r$i"))
        else
            mi[mi.<qx[1]] .= qx[1]
            mi[mi.>qx[2]] .= qx[2]
            mi = (mi .- qx[1]) ./ (qx[2] - qx[1])
            b1, b2 = rvals[i] isa Number ? (rvals[i], rvals[i]) : rvals[i]
            m[i, :] .= (b1 .+ (b2 - b1) .* mi)
        end
    end
    lp
end

function rescale_magnitude(lpo::LocalAnisotropy; kwargs...)
    lp = deepcopy(lpo)
    rescale_magnitude!(lp; kwargs...)
end

"""
    reference_magnitude(localaniso, axis)
    reference_magnitude!(localaniso, axis)

Rescale the magnitudes of `localaniso` by setting `axis` as a unit reference.
`axis` can be :X, :Y or :Z.

## Example

For example, consider a point with magnitude equal to [1.0, 0.2]. After scaling
it with a ExponentialVariogram of range=10m, the range in X would be 10m and in
Y would be 2m at this point. If we first call `reference_magnitude!(localaniso, :Y)`,
the magnitude will turn to [5.0, 1.0], and after scaling to the same previous
variogram, the range in X would be 50m and in Y would be 10m.
"""
function reference_magnitude!(lpar::LocalAnisotropy, ax::Symbol)
    # rescale magnitude according to reference axis
    m = lpar.magnitude
    ref = m[iaxis(ax), :]
    for i = 1:size(m, 1)
        m[i, :] ./= ref
    end
    lpar
end

function reference_magnitude(lpar::LocalAnisotropy, ax::Symbol)
    lp = deepcopy(lpar)
    reference_magnitude!(lp, ax)
end

"""
    adjust_rake(localaniso, azimuths)
    adjust_rake!(localaniso, azimuths)

Adjust the main direction :X along the azimuth informed.
This will change the directions of :X and :Y, and keep the :Z unchanged
Useful when the ellipsoid has the correct :Z (e.g. along body thickness), but the
main direction should align with an specific azimuth onto the main plane.

## Example

Considering three ellipsoids in the `localaniso` variable

```julia
azimuths = [25.0, 30.0, 40.0]
newpars = adjust_rake(localaniso, azimuths)
```
"""
function adjust_rake!(lpar::LocalAnisotropy, az::AbstractVector)

    rot = map(1:nvals(lpar)) do i
        dcm = quat_to_dcm(lpar.rotation[i])
        p1 = Plane(Point(0, 0, 0), Vec(dcm[1, :]), Vec(dcm[2, :]))
        n2 = Vec(sind(az[i] + 90), cosd(az[i] + 90), 0)
        p2 = Plane(Point(0, 0, 0), n2)

        dcm = collect(dcm)
        icross = ustrip.(to(intersection(p1, p2).geom.b))
        dcm[1, :] .= icross
        dcm[2, :] .= cross(dcm[1, :], dcm[3, :])
        det(dcm) < 0 && (dcm = Diagonal([-1, 1, 1]) * dcm)
        dcm_to_quat(DCM(dcm))
    end
    lpar.rotation .= rot
    lpar
end

function adjust_rake(lpar::LocalAnisotropy, az::AbstractVector)
    lp = deepcopy(lpar)
    adjust_rake!(lp, az)
end

"""
    to_vtk(vtkfile, coords, localaniso; dir=:ellips, magnitude=:ranges)

Export local anisotropies `localaniso` at `coords` to VTK format. It export local
ellipsoids by default. If `dir` is set to `:X`, `:Y`, `:Z` or `:XYZ`, local
vectors/arrows are exported on these axes. The magnitude can be expressed as
`:ranges` (default) or `:ratios`. In the case of `:ratios`, `ratio1=range2/range1`
and `ratio2=range3/range1`. Instead of an array of coordinates, `coords` can
also be just the georeferenced spatial object. The output `.vtu` file can be
loaded in Paraview or similars and processed as Glyphs to plot the directions.
If `dir = :ellips`, TensorGlyphs must be used for proper visualization (need to
disable extract eigenvalues).
"""
function to_vtk(
    vtkfile,
    coords::AbstractArray,
    lpars::LocalAnisotropy;
    dir = :ellips,
    magnitude = :ranges,
)
    @assert dir in [:ellips, :XYZ, :X, :Y, :Z] "`dir` must be :ellips, :XYZ :X, :Y or :Z"
    @assert magnitude in [:ratios, :ranges] "`magnitude` must be :ratios or :ranges"

    elp = (dir == :ellips)
    dim = size(coords, 1)
    n = size(coords, 2)
    axs = dir in [:XYZ, :ellips] ? (1:dim) : [iaxis(dir)]
    xyz = dim == 2 ? vcat(coords, zeros(Float64, 1, n)) : coords
    vtx = [MeshCell(VTKCellTypes.VTK_VERTEX, [i]) for i = 1:n]
    ijk = elp ? zeros(Float64, 9, n, 1, 1) : [zeros(Float64, 3, n, 1, 1) for i in axs]
    zx = (dim == 2 && elp) ? minimum(lpars.magnitude) / 10 : 0.0

    # put data into vtk structure format
    for x = 1:n
        if elp
            mag = lpars.magnitude[:, x]
            dim == 2 && push!(mag, zx)
            Λ = Diagonal(SVector{3}(mag))
            P = rotmat(lpars, x)' * Λ
            for v = 1:length(P)
                ijk[v, x, 1, 1] = P[v]
            end
        else
            dcm = rotmat(lpars, x)
            f = dcm[3, 3] < 0 ? -1 : 1
            for d in axs
                ijk[d][1, x, 1, 1] = f * dcm[d, 1]
                ijk[d][2, x, 1, 1] = f * dcm[d, 2]
                ijk[d][3, x, 1, 1] = f * dcm[d, 3]
            end
        end
    end

    # write vtk file
    outfiles = vtk_grid(vtkfile, xyz, (vtx)) do vtk
        if elp
            vtk["ellips"] = ijk
        else
            for d in axs
                vtk["dir$d"] = ijk[d]
            end
        end
        for d = 1:dim
            if magnitude == :ranges
                vtk["range$d"] = lpars.magnitude[d, :]
            elseif d > 1
                vtk["ratio$(d-1)"] = lpars.magnitude[d, :] ./ lpars.magnitude[1, :]
            end
        end
    end
end

function to_vtk(vtkfile, D::SpatialData, lpars::LocalAnisotropy; kwargs...)
    coords = reduce(hcat, [ustrip.(to(centro(D, x))) for x = 1:nvals(D)])
    to_vtk(vtkfile, coords, lpars; kwargs...)
end

toqmat(lp) =
    isnothing(lp) ? nothing : [qmat(rotation(lp, i), magnitude(lp, i)) for i = 1:nvals(lp)]

Base.vcat(lpars::LocalAnisotropy...; kwars...) = reduce(vcat, lpars)

function Base.vcat(lpar1::LocalAnisotropy, lpar2::LocalAnisotropy; kind = :union)
    kind != :union && error("vcat of LocalAnisotropy is restricted to kind=:union ")
    rot = vcat(lpar1.rotation, lpar2.rotation)
    mag = hcat(lpar1.magnitude, lpar2.magnitude)
    LocalAnisotropy(rot, mag)
end
