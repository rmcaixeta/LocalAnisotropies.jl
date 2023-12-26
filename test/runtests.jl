using LocalAnisotropies
using GeoStats
using Test
import LocalAnisotropies: rotmat


@testset "LocalAnisotropies.jl" begin

    # convert data to LocalAnisotropy
    dummy = georef((az=1:10, r1=1:10, r2=1:10), PointSet(reshape(1:20,(2,10))/10))
    pars  = localanisotropies(dummy, [:az], [:r1,:r2], :GSLIB)
	@test round(pars.rotation[1][4],digits=4) ≈ 0.0087

    # convert between different rotation conventions
    angs1 = convertangles([30,30,30], :GSLIB, :Datamine)
    angs2 = convertangles(angs1, :Datamine, :GSLIB)
	rmat  = rotmat([1,1,1], [30,30,30], :GSLIB)
    @test all(rmat .≈ rotmat([1,1,1], angs1, :Datamine))
    @test all(rmat .≈ rotmat([1,1,1], angs2, :GSLIB))
    angs_ = convertangles.(pars.rotation, :GSLIB)
	@test all([x[1] for x in angs_] .≈ 1:10)

    # pars_ = convertpars(pars, convention=:Datamine)
    # exportpars("C:\\Users\\test.csv", pars, :GSLIB)

	# local anisotropies from pointset coordinates
	data  = rmat[1] * vcat(reshape(1:200,(2,100))/10,zeros(1,100))
	pset  = PointSet(data)
	nhood = KNearestSearch(pset, 10)
    gpars = localanisotropies(Geometrical, nhood, simplify=true)
    gpars = localanisotropies(Geometrical, nhood, simplify=false)

    grid2d = (10,10)
    grid3d = (10,10,5)

    for dims in (grid2d, grid3d)
        # reference scenario for tests
		R = length(dims)
		D = if R==2
    		georef((P=[sin(i)+j for i in 1:dims[1], j in 1:dims[2]],))
    	else
    		georef((P=[sin(i)+j+k for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]],))
    	end

        n = round(Int, 0.2*prod(dims))
        S = sample(D, UniformSampling(n, replace=false))
        S = georef(values(S), PointSet(centroid.(domain(S))))
        G = CartesianGrid(dims...)
        searcher = KNearestSearch(G, 10)

        # get local anisotropies
        lpars = localanisotropies(Gradients, D, :P, 3)

        # rescale magnitude and interpolate local anisotropies
        lpars = rescale_magnitude(lpars, r2=(0.2,1.0))
        lpars = reference_magnitude(lpars, :Y)
        lpars = smoothpars(lpars, searcher)

		# pass to samples and calculate expvario
		spars = nnpars(lpars, D, S)
		expvario = localvariography(S, spars, :P, tol=2, maxlag=20, nlags=20, axis=:X)

        # interpolate in a coarser grid
        G_ = if R==2
			CartesianGrid((5,5),(0.5,0.5),(2.0,2.0))
		else
			CartesianGrid((5,5,3),(0.5,0.5,0.5),(2.0,2.0,2.0))
		end
        lpars_ = idwpars(lpars, searcher, G_)

        # Variogram
        γ = NuggetEffect(0.1) + 0.9*ExponentialVariogram(range=60.0)

        # LocalKriging (MW)
        MW = LocalKriging(variogram=γ, localaniso=lpars, method=:MovingWindows)
		s1 = S |> InterpolateNeighbors(G, :P=>MW, maxneighbors=20)

        # LocalKriging (KC)
        KC = LocalKriging(variogram=γ, localaniso=lpars, method=:KernelConvolution)
        s1 = S |> InterpolateNeighbors(G, :P=>KC, maxneighbors=20)

        # Spatial deformation: anisotropic distances
        Sd1, Dd1 = deformspace(S, G, lpars, AnisoDistance, anchors=250)
        s3 = Sd1 |> Interpolate(Dd1, :P=>Kriging(γ))
		x3 = to_3d(s3)

        # Spatial deformation: anisotropic variogram distances
        Sd2, Dd2 = deformspace(S, G, lpars, KernelVariogram, γ, anchors=250)
        s4 = Sd2 |> Interpolate(Dd2, :P=>Kriging(γ))
		x4 = to_3d(s4)

        # Spatial deformation: geodesic anisotropic distances
        LDa = graph(S, G, lpars, AnisoDistance, searcher)
        Sd3, Dd3 = deformspace(LDa, GraphDistance, anchors=250)
        s5 = Sd3 |> Interpolate(Dd3, :P=>Kriging(γ))
		x5 = to_3d(s5)

        # Spatial deformation: geodesic anisotropic variogram distances
        LDv = graph(S, G, lpars, KernelVariogram, γ, searcher)
        Sd4, Dd4 = deformspace(LDv, GraphDistance, anchors=250)
		s6 = Sd4 |> Interpolate(Dd4, :P=>Kriging(γ))
		x6 = to_3d(s6)
    end
end
