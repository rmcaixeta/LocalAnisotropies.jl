
# output conversion main format? a1, a2, a3, r1, r2, r3 (centroid too?)
# extra conventions: MainPlaneNormalAngles, MainPlaneNormal, MainPlaneAngles, MainPlaneVector
# only convert: Dir1Angles, Dir1Vector, ....

# convert angles + ranges into LocalParameters
function LocalParameters(data::GeoData, angles::AbstractVector,
    semiaxes::AbstractVector, convention=:GSLIB)
    semiaxes = Symbol.(semiaxes)
    angles   = Symbol.(angles)
    tab  = values(data)
    cols = string.(propertynames(tab))
    len  = size(tab,1)
    @assert string.(angles)   ⊆ cols "angle column name do not exist"
    @assert string.(semiaxes) ⊆ cols "semiaxis column name do not exist"

    q = Array{Quaternion}(undef,len)
    m = Array{Vector}(undef,len)

    for i in 1:len
        xsemiaxes = [getproperty(tab[i],x) for x in semiaxes] ./ getproperty(tab[i],semiaxes[1])
        xangles   = [getproperty(tab[i],x) for x in angles]

        P, Λ = rotmat(xsemiaxes, xangles, convention; rev=false)
        size(P,1) == 2 && (P=DCM([P[1,1] P[1,2] 0; P[2,1] P[2,2] 0; 0 0 1]))
        q[i] = dcm_to_quat(P)
        m[i] = xsemiaxes
    end

    LocalParameters(q, reduce(hcat,m))
end

# convert angles between different conventions
function convertangles(angles::AbstractArray, convention1::Symbol, convention2::Symbol)
    P, _  = rotmat([1,1,1], angles, convention1)
    rotmat2angles(P, convention2)
end

# convert quaternion to angles in different conventions
function convertangles(quat::Quaternion, convention::Symbol)
    dcm    = quat_to_dcm(quat)
    rotmat2angles(dcm, convention)
end

# convert export local parameters using desired output convention
function exportpars(out, lpars::LocalParameters, convention=:GSLIB)
    table = convertpars(lpars, convention)
    CSV.write(out, table)
end

# convert LocalParameters into angles + ranges ### Simplify it
function convertpars(lpars::LocalParameters, convention=:GSLIB)
    pars = []
    len  = nvals(lpars)
    for i in 1:len
        dcm    = quat_to_dcm(rotation(lpars, i))
        angles = rotmat2angles(dcm, convention)
        ranges = magnitude(lpars, i)
        push!(pars, append!(angles,ranges))
    end
    pars = reduce(hcat, pars)
    r3   = size(pars,1) == 6 ? pars[6,:] : [0 for i in 1:len]
    ((ang1=pars[1,x],ang2=pars[2,x],ang3=pars[3,x],r1=pars[4,x],r2=pars[5,x],r3=r3[x]) for x in 1:len)
end

# reverse transformation of rotation matrix to angles
function rotmat2angles(dcm::AbstractMatrix, convention::Symbol)
  N = size(dcm, 1)
  P = N == 2 ? DCM([P[1,1] P[1,2] 0; P[2,1] P[2,2] 0; 0 0 1]) : dcm

  rule = rules[convention]
  !rule.extrinsic && (P = P')
  preangs = dcm_to_angle(P, rule.order)
  angles  = [preangs.a1, preangs.a2, preangs.a3]

  intr = @. (rule.motion == :CW)  & !rule.extrinsic
  extr = @. (rule.motion == :CCW) & rule.extrinsic
  angles[intr .| extr] *= -1

  !rule.radian && (angles = rad2deg.(angles))
  angles
end
