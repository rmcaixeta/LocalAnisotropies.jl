using Distances
using GeoStats
using LocalAnisotropies
using Test
import LocalAnisotropies: rotmat, anisodistance

@testset "LocalAnisotropies.jl" begin
  d₁ = anisodistance([1.0, 1.0], [0.0])
  d₂ = anisodistance([1.0, 2.0], [0.0])
  @test evaluate(d₁, [1.0, 0.0], [0.0, 0.0]) == evaluate(d₁, [0.0, 1.0], [0.0, 0.0])
  @test evaluate(d₂, [1.0, 0.0], [0.0, 0.0]) != evaluate(d₂, [0.0, 1.0], [0.0, 0.0])

  d₃ = anisodistance([1.0, 0.5, 0.5], [π / 4, 0.0, 0.0])
  @test evaluate(d₃, [1.0, 1.0, 0.0], [0.0, 0.0, 0.0]) ≈ √2
  @test evaluate(d₃, [-1.0, 1.0, 0.0], [0.0, 0.0, 0.0]) ≈ √8

  # intrinsic conventions
  gslib = anisodistance([50.0, 25.0, 5.0], [30.0, -30.0, 30.0], :GSLIB)
  tait = anisodistance([25.0, 50.0, 5.0], [-π / 6, -π / 6, π / 6], :TaitBryanIntr)
  euler = anisodistance([50.0, 25.0, 5.0], [-deg2rad(78), -deg2rad(41), -deg2rad(50)], :EulerIntr)
  lpf = anisodistance([50.0, 25.0, 5.0], [78.0, 41.0, 50.0], :Leapfrog)
  dm = anisodistance([50.0, 25.0, 5.0], [78.0, 41.0, 50.0], :Datamine)

  @test evaluate(gslib, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]) ≈ 0.1325707358356285
  @test evaluate(gslib, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]) ≈ 0.039051248379533283
  @test evaluate(gslib, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]) ≈ 0.15132745950421558

  @test evaluate(gslib, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]) ≈ evaluate(tait, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
  @test evaluate(euler, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]) ≈ evaluate(dm, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
  @test evaluate(euler, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]) ≈ evaluate(lpf, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])
  @test evaluate(euler, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]) - evaluate(gslib, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]) < 10^-3

  # extrinsic conventions
  xtait = anisodistance([50.0, 25.0, 5.0], [π, 0, π / 2], :TaitBryanExtr)
  xeuler = anisodistance([50.0, 25.0, 5.0], [-π / 2, -π / 2, -π / 2], :EulerExtr)

  @test evaluate(xtait, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]) ≈ 0.20
  @test evaluate(xtait, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]) ≈ 0.04
  @test evaluate(xtait, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]) ≈ 0.02

  @test evaluate(xtait, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]) ≈ evaluate(xeuler, [1.0, 0.0, 0.0], [0.0, 0.0, 0.0])
  @test evaluate(xtait, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]) ≈ evaluate(xeuler, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0])
  @test evaluate(xtait, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0]) ≈ evaluate(xeuler, [0.0, 0.0, 1.0], [0.0, 0.0, 0.0])

  # convert data to LocalAnisotropy
  dummy = georef((az=1:10, r1=1:10, r2=1:10), PointSet([(i / 2, (i + 1) / 2) for i in 1:2:19]))
  pars = localanisotropies(dummy, [:az], [:r1, :r2], :GSLIB)

  tabpars_quat = to_table(pars)
  tabpars_angs = to_table(pars, :GSLIB)
  tabpars_quat_xyz = to_table(domain(dummy), pars)
  tabpars_angs_xyz = to_table(domain(dummy), pars, :GSLIB)

  @test length(keys(tabpars_quat)) == 6
  @test length(keys(tabpars_angs)) == 5
  @test length(keys(tabpars_quat_xyz)) == 8
  @test length(keys(tabpars_angs_xyz)) == 7

  repars = localanisotropies(tabpars_quat, [:r1, :r2])
  repars = localanisotropies(tabpars_angs, [:r1, :r2], :GSLIB)

  # convert between different rotation conventions
  angs1 = convertangles([30, 30, 30], :GSLIB, :Datamine)
  angs2 = convertangles(angs1, :Datamine, :GSLIB)
  rmat, _ = rotmat([1, 1, 1], [30, 30, 30], :GSLIB)
  dmat, _ = rotmat([1, 1, 1], angs1, :Datamine)
  gmat, _ = rotmat([1, 1, 1], angs2, :GSLIB)
  @test all(rmat .≈ dmat)
  @test all(rmat .≈ gmat)
  angs_ = convertangles.(pars.rotation, :GSLIB)
  @test all([x[1] for x in angs_] .≈ 1:10)

  # convert data to 3D LocalAnisotropy
  pts3d = rand(3, 10)
  dummy =
    georef((a1=1:10, a2=1:10, a3=1:10, r1=1:10, r2=1:10, r3=1:10), PointSet([tuple(pts3d[:, x]...) for x in 1:10]))
  pars = localanisotropies(dummy, [:a1, :a2, :a3], [:r1, :r2, :r3], :Datamine)
  pars_ = adjust_rake(pars, [45 + i for i in 1:10])
  @test all(round.(rotmat(pars, 1)[3, :], digits=4) .≈ round.(rotmat(pars_, 1)[3, :], digits=4))
  @test !all(round.(rotmat(pars, 1)[1, :], digits=4) .≈ round.(rotmat(pars_, 1)[1, :], digits=4))

  # try some utils
	normals = to_vector(pars, 3)
	azplunge_maindir = to_vector_angle(pars, 1)
	dipdir_dip = to_plane_angle(pars)
	dipvector = to_plane_dipvector(pars)

  tabpars_quat = to_table(pars)
  tabpars_angs = to_table(pars, :Datamine)
  tabpars_quat_xyz = to_table(domain(dummy), pars)
  tabpars_angs_xyz = to_table(domain(dummy), pars, :Datamine)

  @test length(keys(tabpars_quat)) == 7
  @test length(keys(tabpars_angs)) == 6
  @test length(keys(tabpars_quat_xyz)) == 10
  @test length(keys(tabpars_angs_xyz)) == 9

  repars = localanisotropies(tabpars_quat, [:r1, :r2, :r3])
  repars = localanisotropies(tabpars_angs, [:r1, :r2, :r3], :Datamine)

  # local anisotropies from pointset coordinates
  data = rmat[1] * vcat(reshape(1:200, (2, 100)) / 10, zeros(1, 100))
  pset = PointSet([tuple(data[:, x]...) for x in 1:100])
  nhood = KNearestSearch(pset, 10)
  gpars = localanisotropies(Geometric, nhood, simplify=true)
  gpars = localanisotropies(Geometric, nhood, simplify=false)

  grid2d = (10, 10)
  grid3d = (10, 10, 5)

  for dims in (grid2d, grid3d)
    # reference scenario for tests
    R = length(dims)
    D = if R == 2
      georef((P=[sin(i) + j for i in 1:dims[1], j in 1:dims[2]],))
    else
      georef((P=[sin(i) + j + k for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]],))
    end

    n = round(Int, 0.2 * prod(dims))
    S = sample(D, UniformSampling(n, replace=false))
    S = georef(values(S), PointSet(centroid.(domain(S))))
    G = CartesianGrid(dims...)
    searcher = KNearestSearch(G, 10)

    # get local anisotropies
    lpars = localanisotropies(Gradients, D, :P, 3)

    # rescale magnitude and interpolate local anisotropies
    lpars = rescale_magnitude(lpars, r2=(0.2, 1.0))
    lpars = reference_magnitude(lpars, :Y)
    lpars = smoothpars(lpars, searcher)

    # pass to samples and calculate expvario
    spars = nnpars(lpars, D, S)
    expvario = localvariography(S, spars, :P, tol=2, maxlag=20, nlags=20, axis=:X)

    # interpolate in a coarser grid
    G_ = if R == 2
      CartesianGrid((5, 5), (0.5, 0.5), (2.0, 2.0))
    else
      CartesianGrid((5, 5, 3), (0.5, 0.5, 0.5), (2.0, 2.0, 2.0))
    end
    lpars_ = idwpars(lpars, searcher, G_)

    # Variogram
    γ = NuggetEffect(0.1) + 0.9 * ExponentialVariogram(range=60.0)

    # LocalIDW (IDW)
    println("IDW")
    ID = LocalIDW(2.0, lpars)
    S |> LocalInterpolate(G, model=ID, maxneighbors=20)

    # Global
    println("GLOBAL")
    GB = LocalKriging(:Global, lpars, γ)
    S |> LocalInterpolate(G, model=GB, maxneighbors=20)

    # LocalKriging (MW)
    println("MW")
    MW = LocalKriging(:MovingWindows, lpars, γ)
    S |> LocalInterpolate(G, model=MW, maxneighbors=20)

    # LocalKriging (KC)
    println("KC")
    KC = LocalKriging(:KernelConvolution, lpars, γ)
    S |> LocalInterpolate(G, model=KC, maxneighbors=6)

    # Spatial deformation: anisotropic distances
    println("SD")
    Sd1, Dd1 = deformspace(S, G, lpars, AnisoDistance, anchors=250)
    s3 = Sd1 |> Interpolate(Dd1, model=Kriging(γ))
    to_3d(s3)

    # Spatial deformation: anisotropic variogram distances
    Sd2, Dd2 = deformspace(S, G, lpars, KernelVariogram, γ, anchors=250)
    s4 = Sd2 |> Interpolate(Dd2, model=Kriging(γ))
    to_3d(s4)

    # Spatial deformation: geodesic anisotropic distances
    LDa = graph(S, G, lpars, AnisoDistance, searcher)
    Sd3, Dd3 = deformspace(LDa, GraphDistance, anchors=250)
    s5 = Sd3 |> Interpolate(Dd3, model=Kriging(γ))
    to_3d(s5)

    # Spatial deformation: geodesic anisotropic variogram distances
    LDv = graph(S, G, lpars, KernelVariogram, γ, searcher)
    Sd4, Dd4 = deformspace(LDv, GraphDistance, anchors=250)
    s6 = Sd4 |> Interpolate(Dd4, model=Kriging(γ))
    to_3d(s6)

    # Local sequential gaussian simulation
    println("SGS")
    local_sgs = LocalSGS(lpars, maxneighbors=10)
    rand(GaussianProcess(γ), G, 2; method=local_sgs, data=S |> Quantile())

    local_sgs = LocalSGS(:KernelConvolution, lpars, maxneighbors=10)
    rand(GaussianProcess(γ), G, 2; method=local_sgs, data=S |> Quantile())
  end
end
