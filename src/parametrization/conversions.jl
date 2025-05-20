# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    localanisotropies(data, rotation, ranges, convention=nothing)
    localanisotropies(data, ranges, convention) # assumes angles [:ang1,:ang2,:ang3]
    localanisotropies(data, ranges) # assumes quaternions [:q0,:q1,:q2,:q3]

Import a `LocalAnisotropy` object from tabular data. The `data` must follow
`Tables.jl` structure. It must have the rotation info and ranges as
properties. These property names are informed as vectors in `rotation` and
`ranges`for the conversion. Check out the available rotation conventions at
[`RotationRule`](@ref) docstring.

## Example

```julia
localanisotropies(data, [:rot1], [:range1, :range2], :EulerIntr) # 2D
localanisotropies(data, [:rot1, :rot2, :rot3], [:range1, :range2, :range3], :GSLIB) # 3D
localanisotropies(data, [:q0,:q1,:q2,:q3], [:range1, :range2, :range3]) # quaternion
```
"""
localanisotropies(data, ranges::AbstractVector) = localanisotropies(data, [:q0, :q1, :q2, :q3], ranges, nothing)

localanisotropies(data, ranges::AbstractVector, convention::RotConvention) =
  localanisotropies(data, [:ang1, :ang2, :ang3], ranges, convention)

function localanisotropies(
  data,
  rotation::AbstractVector,
  ranges::AbstractVector,
  convention::Union{RotConvention,Nothing}=nothing
)
  ranges = Symbol.(ranges)
  rotation = Symbol.(rotation)
  vals = data isa NamedTuple ? data : values(data)
  tab = Tables.columns(vals)
  cols = string.(Tables.columnnames(tab))
  len = nrow(data)

  @assert string.(rotation) ⊆ cols "angle column name do not exist"
  @assert string.(ranges) ⊆ cols "range column name do not exist"
  dim = length(ranges)
  transf = !isnothing(convention)
  isquat = length(rotation) == 4

  q = Array{Quaternion}(undef, len)
  m = Array{Vector}(undef, len)

  if transf
    rule = convention isa Symbol ? rules[convention] : convention
    if rule.main == :y
      ranges = ranges[reverse(1:dim, 1, 2)]
      rule = RotationRule(rule.order, rule.motion, rule.radian, :x, rule.extrinsic)
    end
  end

  for i in 1:len
    vranges = [Tables.getcolumn(tab, x)[i] for x in ranges]
    xranges = vranges ./ maximum(vranges)
    xrot = [Tables.getcolumn(tab, x)[i] for x in rotation]

    quat = if isquat
      Quaternion(xrot)
    else
      P, Λ = rotmat(xranges, xrot, rule)
      size(P, 1) == 2 && (P = DCM([P[1, 1] P[1, 2] 0; P[2, 1] P[2, 2] 0; 0 0 1]))
      dcm_to_quat(P)
    end

    q[i] = quat
    m[i] = xranges
  end

  LocalAnisotropy(q, reduce(hcat, m))
end

function localanisotropies(vectors::AbstractArray; dipvector=false)
  # only adapted to 3d
  conversion_func = dipvector ? dipvector_to_quaternion : normal_to_quaternion
  q = if size(vectors, 1) == 3
    size(vectors, 2) == 3 && println("Warning: 3x3 table might return strange results if vectors are not in column-major order.")
    conversion_func.(eachcol(vectors))
  else
    conversion_func.(eachrow(vectors))
  end

  LocalAnisotropy(q, ones(3, length(q)))
end

"""
    convertangles(angles, convention1, convention2) # angles to new angles
    convertangles(angles, convention1) # angles to quaternion

Convert angles between different convention. Check out the available rotation
conventions at [`RotationRule`](@ref) docstring.

## Example

Converting 3D rotation [30, 30 30] from GSLIB convention to Datamine convention

```julia
new_angs = convertangles([30,30,30], :GSLIB, :Datamine)
```
"""
function convertangles(angles::AbstractVector, c1::RotConvention, c2::RotConvention)
  P, _ = rotmat([1, 1, 1], angles, c1)
  rotmat2angles(P, c2)
end

function convertangles(angles::AbstractVector, c1::RotConvention)
  P, _ = rotmat([1, 1, 1], angles, c1)
  size(P, 1) == 2 && (P = DCM([P[1, 1] P[1, 2] 0; P[2, 1] P[2, 2] 0; 0 0 1]))
  dcm_to_quat(P)
end

"""
    convertangles(quaternion, convention)

Convert a `ReferenceFrameRotations.Quaternion` to angles in a given convention.
Check out the available rotation conventions at [`RotationRule`](@ref) docstring.
"""
function convertangles(quat::Quaternion, convention::RotConvention)
  dcm = quat_to_dcm(quat)
  rotmat2angles(dcm, convention)
end

"""
    to_table(localaniso)
    to_table(domain, localaniso)
    to_table(localaniso, convention)
    to_table(domain, localaniso, convention)

Transform local anisotropies to a Table file. The angles are written in the form of
the given convention (or as quaternions if not informed). Check out the available
rotation conventions at [`RotationRule`](@ref) docstring. Can be materialized as
DataFrame, form example, or kept as NamedTuple format.

## Example

```julia
using DataFrames, CSV
to_table(localaniso, :GSLIB) |> DataFrame

CSV.write("out.csv",to_table(localaniso))
```
"""
# convert LocalGeoData into quaternions + ranges
to_table(obj::SpatialData, lpars::LocalAnisotropy) = to_table(LocalGeoData(obj, lpars))

function to_table(lpars::LocalGeoData)
  qs = [Symbol("q$(i-1)") => [x[i] for x in rotation(lpars)] for i in 1:4]
  mag = [Symbol("r$i") => magnitude(lpars)[i, :] for i in 1:ndims(lpars)]
  cnames = [:x, :y, :z]
  cvals = coords_(obj(lpars))
  cvals = [cnames[i] => cvals[i, :] for i in 1:ndims(lpars)]
  (; cvals..., qs..., mag...)
end

# convert LocalAnisotropy into quaternions + ranges
function to_table(lpars::LocalAnisotropy)
  qs = [Symbol("q$(i-1)") => [x[i] for x in rotation(lpars)] for i in 1:4]
  mag = [Symbol("r$i") => magnitude(lpars)[i, :] for i in 1:ndims(lpars)]
  (; qs..., mag...)
end

# convert LocalAnisotropy into angles + ranges
function to_table(lpars::LocalAnisotropy, convention::RotConvention)
  pars = mapreduce(hcat, 1:nvals(lpars)) do i
    dcm = rotmat(lpars, i)
    angles = rotmat2angles(dcm, convention)
    ranges = magnitude(lpars, i)
    vcat(angles, ranges)
  end
  cols = [:ang1, :ang2, :ang3, :r1, :r2, :r3]
  pairs = [cols[i] => pars[i, :] for i in 1:size(pars, 1)]
  (; pairs...)
end

function to_table(lgeo::LocalGeoData, convention::RotConvention)
  lpars = lgeo.localaniso
  pars = mapreduce(hcat, 1:nvals(lpars)) do i
    dcm = rotmat(lpars, i)
    angles = rotmat2angles(dcm, convention)
    ranges = magnitude(lpars, i)
    vcat(angles, ranges)
  end
  cols = [:ang1, :ang2, :ang3, :r1, :r2, :r3]
  pairs = [cols[i] => pars[i, :] for i in 1:size(pars, 1)]

  cnames = [:x, :y, :z]
  cvals = coords_(obj(lgeo))
  cvals = [cnames[i] => cvals[i, :] for i in 1:ndims(lgeo)]
  (; cvals..., pairs...)
end

to_table(obj::SpatialData, lpars::LocalAnisotropy, convention::RotConvention) =
  to_table(LocalGeoData(obj, lpars), convention)

# deprecated
function convertpars(args...)
  @warn "convertpars deprecated; prefer using to_table function"
  to_table(args...)
end

# reverse transformation of rotation matrix to angles
function rotmat2angles(dcm::AbstractMatrix, convention::RotConvention)
  N = size(dcm, 1)
  P = N == 2 ? DCM([P[1, 1] P[1, 2] 0; P[2, 1] P[2, 2] 0; 0 0 1]) : dcm

  rule = convention isa Symbol ? rules[convention] : convention
  rule.extrinsic && (P = P')
  preangs = dcm_to_angle(DCM(P), rule.order)
  angles = [preangs.a1, preangs.a2, preangs.a3]

  intr = @. (rule.motion == :CW) & !rule.extrinsic
  extr = @. (rule.motion == :CCW) & rule.extrinsic
  angles[intr .| extr] *= -1

  !rule.radian && (angles = rad2deg.(angles))
  angles
end
