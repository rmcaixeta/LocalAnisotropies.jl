# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    localanisotropies(data, angles, ranges, convention=:GSLIB)

Import a `LocalAnisotropy` object from an outer source. The `data` is a
georeferenced object. It must have the rotation angles and ranges/ratios as
properties. These property names are informed as vectors in `angles` and
`ranges`for the conversion. Check out the available rotation conventions at
[`RotationRule`](@ref) docstring.

## Example

```julia
localanisotropies(data, [:rot1], [:range1, :range2], convention=:EulerIntr) # 2D
localanisotropies(data, [:rot1, :rot2, :rot3], [:range1, :range2, :range3]) # 3D
```
"""
function localanisotropies(data::SpatialData, angles::AbstractVector,
    ranges::AbstractVector, convention=:GSLIB)
    ranges = Symbol.(ranges)
    angles = Symbol.(angles)
    tab  = Tables.columns(values(data))
    cols = string.(Tables.columnnames(tab))
    len  = nelements(data)
    @assert string.(angles) ⊆ cols "angle column name do not exist"
    @assert string.(ranges) ⊆ cols "range column name do not exist"

    q = Array{Quaternion}(undef,len)
    m = Array{Vector}(undef,len)

    for i in 1:len
        xranges = [Tables.getcolumn(tab,x)[i] for x in ranges] ./ Tables.getcolumn(tab,ranges[1])[i]
        xangles = [Tables.getcolumn(tab,x)[i] for x in angles]

        P, Λ = rotmat(xranges, xangles, convention; rev=false)
        size(P,1) == 2 && (P=DCM([P[1,1] P[1,2] 0; P[2,1] P[2,2] 0; 0 0 1]))
        q[i] = dcm_to_quat(P)
        m[i] = xranges
    end

    LocalAnisotropy(q, reduce(hcat,m))
end

#
"""
    convertangles(angles, convention1, convention2)

Convert angles between different convention. Check out the available rotation
conventions at [`RotationRule`](@ref) docstring.

## Example

Converting 3D rotation [30, 30 30] from GSLIB convention to Datamine convention

```julia
new_angs = convertangles([30,30,30], :GSLIB, :Datamine)
```
"""
function convertangles(angles::AbstractVector, convention1, convention2)
    P, _  = rotmat([1,1,1], angles, convention1)
    rotmat2angles(P, convention2)
end

"""
    convertangles(quaternion, convention)

Convert a `ReferenceFrameRotations.Quaternion` to angles in a given convention.
Check out the available rotation conventions at [`RotationRule`](@ref) docstring.
"""
function convertangles(quat::Quaternion, convention)
    dcm    = quat_to_dcm(quat)
    rotmat2angles(dcm, convention)
end

"""
    exportpars(filename, localaniso, convention)

Export local anisotropies to a CSV table. The angles are exported in the form of
the given convention. Check out the available rotation conventions at
[`RotationRule`](@ref) docstring.

## Example

```julia
exportpars("path/localaniso.csv", localaniso, :GSLIB)
```
"""
function exportpars(out, lpars::LocalAnisotropy, convention=:GSLIB)
    table = convertpars(lpars, convention)
    CSV.write(out, table)
end

# convert LocalAnisotropy into angles + ranges
function convertpars(lpars::LocalAnisotropy, convention=:GSLIB)
    pars = []
    len  = nvals(lpars)
    for i in 1:len
        dcm    = rotmat(lpars, i)
        angles = rotmat2angles(dcm, convention)
        ranges = magnitude(lpars, i)
        push!(pars, append!(angles,ranges))
    end
    pars = reduce(hcat, pars)
    r3   = size(pars,1) == 6 ? pars[6,:] : [0 for i in 1:len]
    ((ang1=pars[1,x],ang2=pars[2,x],ang3=pars[3,x],r1=pars[4,x],r2=pars[5,x],r3=r3[x]) for x in 1:len)
end

# reverse transformation of rotation matrix to angles
function rotmat2angles(dcm::AbstractMatrix, convention)
  N = size(dcm, 1)
  P = N == 2 ? DCM([P[1,1] P[1,2] 0; P[2,1] P[2,2] 0; 0 0 1]) : dcm

  rule = convention isa Symbol ? rules[convention] : convention
  !rule.extrinsic && (P = P')
  preangs = dcm_to_angle(P, rule.order)
  angles  = [preangs.a1, preangs.a2, preangs.a3]

  intr = @. (rule.motion == :CW)  & !rule.extrinsic
  extr = @. (rule.motion == :CCW) & rule.extrinsic
  angles[intr .| extr] *= -1

  !rule.radian && (angles = rad2deg.(angles))
  angles
end
