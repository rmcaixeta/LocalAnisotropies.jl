# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from an older version of GeoStatsBase.jl
# ------------------------------------------------------------------

"""
    RotationRule(order, motion, radian, main, extrinsic)

Creates a rule for the ellipsoidal rotation

* order     - sequence of three axes by which the rotations are made; e.g. `:ZXZ`
* motion    - inform for each of the three rotations if it is clockwise (`:CW`) or
              counterclockwise (`:CCW`). Right-hand rule: motion defined looking
              towards the negative direction of the axis
* radian    - `true` if the input angles are in radians or `false` if in degrees
* main      - inform if the main semiaxis is `:x` or `:y`
* extrinsic - `true` if rotation is extrinsic or `false` if it is intrinsic

## Example

Rotation rule to reproduce GSLIB rotation

```julia
rule = RotationRule(:ZXY,[:CW,:CCW,:CCW],false,:y,false)
```

## Conventions

Some rotation conventions built-in:

- `:TaitBryanExtr` => Extrinsic right-handed rotation by the ZXY axes
- `:TaitBryanIntr` => Intrinsic right-handed rotation by the ZXY axes
- `:EulerExtr`     => Extrinsic right-handed rotation by the ZXZ axes
- `:EulerIntr`     => Intrinsic right-handed rotation by the ZXZ axes
- `:GSLIB`         => GSLIB software rotation convention
- `:Leapfrog`      => Leapfrog software rotation convention
- `:Datamine`      => Datamine software rotation convention, axes order ZXZ
- `:Datamine313`   => Datamine software rotation convention, axes order ZXZ
- `:Datamine321`   => Datamine software rotation convention, axes order ZYX
- `:Datamine312`   => Datamine software rotation convention, axes order ZXY
- `:Datamine323`   => Datamine software rotation convention, axes order ZYZ
"""
struct RotationRule
  order::Symbol
  motion::Vector{Symbol}
  radian::Bool
  main::Symbol
  extrinsic::Bool
end

# list of available rotation rules
rules = Dict(
  :TaitBryanExtr => RotationRule(:ZXY,[:CCW,:CCW,:CCW],true,:x,true),
  :TaitBryanIntr => RotationRule(:ZXY,[:CCW,:CCW,:CCW],true,:x,false),
  :EulerExtr     => RotationRule(:ZXZ,[:CCW,:CCW,:CCW],true,:x,true),
  :EulerIntr     => RotationRule(:ZXZ,[:CCW,:CCW,:CCW],true,:x,false),
  :GSLIB         => RotationRule(:ZXY,[:CW,:CCW,:CCW],false,:y,false),
  :Leapfrog      => RotationRule(:ZXZ,[:CW,:CW,:CW],false,:x,false),
  :Datamine      => RotationRule(:ZXZ,[:CW,:CW,:CW],false,:x,false),
  :Datamine313   => RotationRule(:ZXZ,[:CW,:CW,:CW],false,:x,false),
  :Datamine321   => RotationRule(:ZYX,[:CW,:CW,:CW],false,:x,false),
  :Datamine312   => RotationRule(:ZXY,[:CW,:CW,:CW],false,:x,false),
  :Datamine323   => RotationRule(:ZYZ,[:CW,:CW,:CW],false,:x,false)
)

RotConvention = Union{Symbol,RotationRule}

function rotmat(semiaxes::AbstractVector, angles::AbstractVector,
                convention=:TaitBryanExtr; rev=false)
  N = length(semiaxes)
  @assert all(semiaxes .> 0) "semiaxes must be positive"
  @assert N ∈ [2,3] "dimension must be either 2 or 3"

  rule = convention isa Symbol ? rules[convention] : convention

  # invert x and y if necessary
  if rule.main == :y
     semiaxes[1], semiaxes[2] = semiaxes[2], semiaxes[1]
  end

  # scaling matrix
  Λ = Diagonal(SVector{N}(one(eltype(semiaxes))./semiaxes.^2))

  # convert to radian and invert sign if necessary
  !rule.radian && (angles = deg2rad.(angles))
  _0 = zero(eltype(angles))
  N == 2 && (angles = [angles[1], _0, _0])
  intr = @. (rule.motion == :CW)  & !rule.extrinsic
  extr = @. (rule.motion == :CCW) & rule.extrinsic
  angles[intr .| extr] *= -1

  # for reverse extrinsic transformation
  rule.extrinsic && rev && (angles *= -1)

  # rotation matrix
  P = angle_to_dcm(angles..., rule.order)[SOneTo(N),SOneTo(N)]
  !rule.extrinsic && (P = P')
  P, Λ
end
