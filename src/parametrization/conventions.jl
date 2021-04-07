# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from an older version of GeoStatsBase.jl
# ------------------------------------------------------------------

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
  :Datamine321   => RotationRule(:ZYZ,[:CW,:CW,:CW],false,:x,false),
  :Datamine312   => RotationRule(:ZXY,[:CW,:CW,:CW],false,:x,false),
  :Datamine323   => RotationRule(:ZYZ,[:CW,:CW,:CW],false,:x,false)
)

function rotmat(semiaxes::AbstractVector, angles::AbstractVector,
                convention::Symbol=:TaitBryanExtr; rev=false)
  N = length(semiaxes)
  @assert all(semiaxes .> 0) "semiaxes must be positive"
  @assert N ∈ [2,3] "dimension must be either 2 or 3"

  rule = rules[convention]

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
