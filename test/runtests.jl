using LocalAnisotropies
using GeoStats
using Test
import LocalAnisotropies: rotmat

@testset "LocalAnisotropies.jl" begin
    # convert data to LocalParameters
    dummy = georef((az=1:10, r1=1:10, r2=1:10), PointSet(rand(2,10)))
    pars  = LocalParameters(dummy, [:az], [:r1,:r2], :GSLIB)

    # convert between different rotation conventions
    angs1 = convertangles([30,30,30], :GSLIB, :Datamine)
    angs2 = convertangles(angs1, :Datamine, :GSLIB)
    @test all(rotmat([1,1,1], [30,30,30], :GSLIB) .≈ rotmat([1,1,1], angs1, :Datamine))
    @test all(rotmat([1,1,1], [30,30,30], :GSLIB) .≈ rotmat([1,1,1], angs2, :GSLIB))
    angs_ = convertangles.(pars.rotation, :GSLIB)

    # pars_ = convertpars(pars, convention=:Datamine)
    # exportpars("C:\\Users\\test.csv", pars, :GSLIB)

    grid2d = (10,10)
    grid3d = (10,10,5)

    for dims in (grid2d, grid3d)
        # reference scenario for tests
        D = if length(dims)==2
    		georef((P=[sin(i)+j for i in 1:dims[1], j in 1:dims[2]],))
    	else
    		georef((P=[sin(i)+j+k for i in 1:dims[1], j in 1:dims[2], k in 1:dims[3]],))
    	end

        n = round(Int, 0.2*prod(dims))
        S = sample(D, n, replace=false)
        G = RegularGrid(dims...)
        searcher = KNearestSearch(G, 10)

        # get local parameters
        lpars = localparameters(Gradients(), D, :P, 3)

        # rescale magnitude and interpolate local parameters
        lpars = rescale_magnitude(lpars, (0.2,1.0))
        lpars = smoothpars(lpars, searcher)

        # interpolate in a coarser grid
        G_ = RegularGrid((5,5),(0.5,0.5),(2.0,2.0))
        lpars_ = IDWpars(lpars, searcher, G_)

        # Estimation problem
        P = EstimationProblem(S, G, :P)
        γ = NuggetEffect(0.1) + 0.9*ExponentialVariogram(range=60.0)

        # LocalKriging (MW)
        MW = LocalKriging(:P => (variogram=(:X=>γ), localpars=lpars, method=:MovingWindows))
        s1 = solve(P, MW)

        # LocalKriging (KC)
        KC = LocalKriging(:P => (variogram=(:X=>γ), localpars=lpars, method=:KernelConvolution))
        s2 = solve(P, KC)

        # Spatial deformation: anisotropic distances
        Sd1, Dd1 = deformspace(S, G, lpars, LocalAnisotropy(), anchors=1500)
        Pd1 = EstimationProblem(Sd1, Dd1, :P)
        s3 = solve(Pd1, Kriging(:P => (variogram=γ,)))

        # Spatial deformation: anisotropic variogram distances
        Sd2, Dd2 = deformspace(S, G, lpars, LocalVariogram(), γ, anchors=1500)
        Pd2 = EstimationProblem(Sd2, Dd2, :P)
        s4 = solve(Pd2, Kriging(:P => (variogram=γ,)))

        # Spatial deformation: geodesic anisotropic distances
        LDa = addgraph(S, G, lpars, LocalAnisotropy(), searcher)
        Sd3, Dd3 = deformspace(LDa, GraphDistance(), anchors=1500)
        Pd3 = EstimationProblem(Sd3, Dd3, :P)
        s5 = solve(Pd3, Kriging(:P => (variogram=γ,)))

        # Spatial deformation: geodesic anisotropic variogram distances
        LDv = addgraph(S, G, lpars, LocalVariogram(), γ, searcher)
        Sd4, Dd4 = deformspace(LDv, GraphDistance(), anchors=1500)
        Pd4 = EstimationProblem(Sd4, Dd4, :P)
        s6 = solve(Pd4, Kriging(:P => (variogram=γ,)))
    end
end
