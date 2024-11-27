### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 714b6cb0-a8ff-11ef-1db8-4743127b2ab7
begin
    import Pkg
    Pkg.activate()
    # load libraries for the example
    using LocalAnisotropies
    using GeoStats
    using Random
    import CairoMakie as Mke
    Random.seed!(1234)

    Mke.set_theme!(size = (500, 500))

    # create a reference scenario for tests
    D = georef((P = [25 - abs(0.2 * i^2 - j) for i = -10:9, j = 1:20],))
    S = sample(D, UniformSampling(80, replace = false))
    S = georef(values(S), PointSet(centroid.(domain(S))))
    G = CartesianGrid(20, 20)

    # plot reference scenario and samples extracted for further estimations
    fig = Mke.Figure(size = (700, 350))
    Mke.plot(fig[1, 1], D.geometry, color = D.P)
    Mke.plot(fig[1, 2], G, showsegments = true, color = :white)
    Mke.plot!(fig[1, 2], S.geometry, color = S.P, pointsize = 10)
    Mke.current_figure()
end

# ╔═╡ f38b20d5-d356-430a-b906-55abbda8a940
begin
    # get local anisotropies using gradients of a 8x8 radius window
    rawlpars = localanisotropies(Gradients, D, :P, 8)
    Mke.plot(D.geometry, color = D.P, alpha = 0.6)
    Mke.plot!(rawlpars, D.geometry)
    Mke.current_figure()
end

# ╔═╡ 52deaf05-3549-49b8-bfe6-0f7c18df44eb
begin
    # average 10 nearest local anisotropies
    searcher = KNearestSearch(G, 10)
    lpars = smoothpars(rawlpars, searcher)

    # rescale magnitude of range2 to between 0.25 and 1.0
    rescale_magnitude!(lpars, r2 = (0.25, 1.0))

    # set possible reference axes to X and Y; ranges will be fixed in those directions
    lparsx = reference_magnitude(lpars, :X)
    lparsy = reference_magnitude(lpars, :Y)

    Mke.plot(D.geometry, color = D.P, alpha = 0.6)
    Mke.plot!(lpars, D.geometry)
    Mke.current_figure()
end

# ╔═╡ 7aa17393-80f2-45e3-ab3d-7dfd4adc2ed0
begin
    # for custom visualizations, it's possible to export it to VTK
    #to_vtk("ellipses", D, lpars)
    # below the file "ellipses.vtu" loaded in Paraview using TensorGlyph (disable extract eigenvalues)
end

# ╔═╡ 534bbe60-cc91-45f4-9e90-7b0d2e576d59
begin
    # pass local anisotropies to samples
    spars = nnpars(lpars, D, S)

    # do an unconventional variography along local X axis (same can be done for Y)
    expvario = localvariography(S, spars, :P, tol = 2, maxlag = 20, nlags = 20, axis = :X)
    Mke.plot(expvario)
    γx = ExponentialVariogram(sill = 32.0, range = 40.0)
    γy = GaussianVariogram(sill = 32.0, range = 8.0)
    Mke.plot!(γx)
    Mke.current_figure()
end

# ╔═╡ 5859126d-a870-49c1-b74f-7eb50e29f78b
begin
    # local inverse distance weighting using exponent 3
    ID = LocalIDW(3, lparsx)
    s1 = S |> LocalInterpolate(G, :P => ID, maxneighbors = 20)
    Mke.plot(s1.geometry, color = s1.P)
end

# ╔═╡ 77869159-c023-4097-8f2f-7af7739f6796
begin
    # kriging using moving windows method
    MW = LocalKriging(:MovingWindows, lparsx, γx)
    s2 = S |> LocalInterpolate(G, :P => MW, maxneighbors = 20)
    Mke.plot(s2.geometry, color = s2.P)
end

# ╔═╡ f8bc4655-f11e-4af9-ae45-d8fdabd99cd2
begin
    # kriging using kernel convolution method (smaller search; unstable with many local variations)
    KC = LocalKriging(:KernelConvolution, lparsy, γy)
    s3 = S |> LocalInterpolate(G, :P => KC, maxneighbors = 6)
    Mke.plot(s3.geometry, color = s3.P)
end

# ╔═╡ 3f48ca1f-8d73-41b4-8721-d19ba881e670
begin
    # deform space using kernel variogram as dissimilarity input
    Sd1, Dd1 = deformspace(S, G, lparsx, KernelVariogram, γx, anchors = 1500)
    γ1 = GaussianVariogram(sill = 21.3, range = 22.5)
    s4 = Sd1 |> Interpolate(Dd1, :P => Kriging(γ1))

    # plot
    fig4 = Mke.Figure(size = (700, 350))
    Mke.plot(fig4[1, 1], to_3d(s4).geometry, color = s4.P)
    Mke.plot(fig4[1, 2], G, color = s4.P)
    Mke.current_figure()
end

# ╔═╡ 85ba6341-92d3-43dd-a9a1-a1e72a60b0d4
begin
    # deform space based on a graph built with average anisotropic distances
    # of the ten nearest data; do variography in multidimensional space
    LDa = graph(S, G, lparsx, AnisoDistance, searcher)
    Sd2, Dd2 = deformspace(LDa, GraphDistance, anchors = 1500)
    γ2 = GaussianVariogram(sill = 22.0, range = 22.0)

    # traditional kriging in the new multidimensional space
    s5 = Sd2 |> Interpolate(Dd2, :P => Kriging(γ2))

    # plot
    fig5 = Mke.Figure(size = (700, 350))
    Mke.plot(fig5[1, 1], to_3d(s5).geometry, color = s5.P)
    Mke.plot(fig5[1, 2], G, color = s5.P)
    Mke.current_figure()
end

# ╔═╡ 5f0f00c9-05cc-4185-aad9-4bb0eab2c1ce
begin
    # deform space based on a graph built with kernel variogram of the ten
    # nearest data; do variography in multidimensional space
    LDv = graph(S, G, lparsy, KernelVariogram, γy, searcher)
    Sd3, Dd3 = deformspace(LDv, GraphDistance, anchors = 1500)
    γ3 = NuggetEffect(1.0) + GaussianVariogram(sill = 21.0, range = 22.0)

    # traditional kriging in the new multidimensional space
    s6 = Sd3 |> Interpolate(Dd3, :P => Kriging(γ3))

    # plot
    fig6 = Mke.Figure(size = (700, 350))
    Mke.plot(fig6[1, 1], to_3d(s6).geometry, color = s6.P)
    Mke.plot(fig6[1, 2], G, color = s6.P)
    Mke.current_figure()
end

# ╔═╡ e935064e-5e18-46ae-919d-83c8ac733f78
begin
    γomni = GaussianVariogram(sill = 32.0, range = 11.0)
    OK = Kriging(γomni)
    s0 = S |> InterpolateNeighbors(G, :P => OK, maxneighbors = 20)
    Mke.plot(s0.geometry, color = s0.P)
end

# ╔═╡ 915af65f-287a-4a11-9bb6-2755bc571a5a
begin
    # comparison of the different estimates
    mse(a, b) = sum((a .- b) .^ 2) / length(b)
    solvers = ["OK", "LIDW", "MW", "KC", "SD1", "SD2", "SD3"]
    errors =
        [mse(getproperty(x, :P), getproperty(D, :P)) for x in [s0, s1, s2, s3, s4, s5, s6]]
    Mke.barplot(
        1:7,
        errors,
        axis = (
            xticks = (1:7, solvers),
            ylabel = "Mean squared error",
            xlabel = "Estimation method",
        ),
    )
end

# ╔═╡ 807f31dc-3e3f-4d2c-ab4a-4eb150f9599b
begin
    # import external local anisotropies in GSLIB convention
    dummy = georef((az = 1:10, r1 = 1:10, r2 = 1:10))
    pars = localanisotropies(dummy, [:az], [:r1, :r2], :GSLIB)

    # interpolate local anisotropies into a coarser grid
    G_ = CartesianGrid((10, 10), (0.5, 0.5), (2.0, 2.0))
    lpars_ = idwpars(lpars, searcher, G_, power = 2.0)

    # convert between different rotation conventions
    angs1 = convertangles([30, 30, 30], :GSLIB, :Datamine)
    angs2 = convertangles.(pars.rotation, :GSLIB)
end

# ╔═╡ 8bb14d5a-1096-4fa8-a505-23b031bb0611
begin
    # normal score transformation
    ns = Quantile()
    S_ns, ref_S_ns = apply(ns, S)
    Sd1_ns, ref_Sd1_ns = apply(ns, Sd1)

    # normal scores unconventional variography along local X axis
    expvario_ns =
        localvariography(S_ns, spars, :P, tol = 2, maxlag = 20, nlags = 20, axis = :X)
    γx_ns = SphericalVariogram(sill = 1.0, range = 15.0)

    # normal scores omni varios
    expomni1_ns = EmpiricalVariogram(S_ns, :P, maxlag = 20, nlags = 20)
    γomni1_ns = SphericalVariogram(sill = 1.0, range = 11.0)
    expomni2_ns = EmpiricalVariogram(Sd1_ns, :P, maxlag = 20, nlags = 20)
    γomni2_ns = ExponentialVariogram(sill = 1.0, range = 50.0)

    figv = Mke.Figure(size = (700, 235))
    Mke.plot(figv[1, 1], expvario_ns)
    Mke.plot!(figv[1, 1], γx_ns, color = :red)
    Mke.plot(figv[1, 2], expomni1_ns)
    Mke.plot!(figv[1, 2], γomni1_ns, color = :red)
    Mke.plot(figv[1, 3], expomni2_ns)
    Mke.plot!(figv[1, 3], γomni2_ns, maxlag = 20, color = :red)
    Mke.current_figure()
end

# ╔═╡ a9dd3f32-92df-4ddf-9f35-6a8f8b73b5be
begin
    # standard sequential gaussian simulation with isotropic variogram
    sgs = SEQMethod(maxneighbors = 25)
    sims1 = rand(GaussianProcess(γomni1_ns), G, S_ns, 100, sgs)
    med1 = quantile(sims1, 0.5)
    sims1_bt = [revert(ns, x, ref_S_ns) for x in (sims1[1], med1)]
    fig7 = Mke.Figure(size = (700, 350))
    Mke.plot(fig7[1, 1], G, color = sims1_bt[1].P, colormap = :jet)
    Mke.plot(fig7[1, 2], G, color = sims1_bt[2].P, colormap = :jet)
    Mke.current_figure()
end

# ╔═╡ e922cb53-d0de-44b9-bb79-c20dfac4263d
begin
    # sequential gaussian simulation with local anisotropies (MW method)
    # not theoretically very correct, but can give good results
    local_sgs = LocalSGS(localaniso = lparsx, maxneighbors = 25)
    sims2 = rand(GaussianProcess(γx_ns), G, S_ns, 100, local_sgs)
    med2 = quantile(sims2, 0.5)
    sims2_bt = [revert(ns, x, ref_S_ns) for x in (sims2[1], med2)]
    fig8 = Mke.Figure(size = (700, 350))
    Mke.plot(fig8[1, 1], G, color = sims2_bt[1].P, colormap = :jet)
    Mke.plot(fig8[1, 2], G, color = sims2_bt[2].P, colormap = :jet)
    Mke.current_figure()
end

# ╔═╡ 48d03eee-4972-443f-a0c9-9fdaafd73c73
begin
    # sequential gaussian simulation after spatial deformation
    # standard simulation after local anisotropies is "removed"
    sims3 = rand(GaussianProcess(γomni2_ns), Dd1, Sd1_ns, 100, sgs)
    med3 = quantile(sims3, 0.5)
    sims3_bt = [revert(ns, x, ref_Sd1_ns) for x in (sims3[1], med3)]
    fig9 = Mke.Figure(size = (700, 350))
    Mke.plot(fig9[1, 1], G, color = sims3_bt[1].P, colormap = :jet)
    Mke.plot(fig9[1, 2], G, color = sims3_bt[2].P, colormap = :jet)
    Mke.current_figure()
end

# ╔═╡ 21ea6b12-e2cb-46df-aed9-0ec6cd2e0cf8
begin
    # comparison of the different simulation methods
    sim_solvers = ["SGS", "SGS_MW", "SGS_SD1"]
    sim_errors = [
        mse(getproperty(x, :P), getproperty(D, :P)) for
        x in [sims1_bt[2], sims2_bt[2], sims3_bt[2]]
    ]
    fig10 = Mke.Figure(size = (700, 350))
    Mke.barplot(
        fig10[1, 1],
        1:3,
        sim_errors,
        axis = (
            xticks = (1:3, sim_solvers),
            ylabel = "Mean squared error",
            xlabel = "Simulation method (median)",
        ),
    )
    Mke.current_figure()
end

# ╔═╡ Cell order:
# ╠═714b6cb0-a8ff-11ef-1db8-4743127b2ab7
# ╠═f38b20d5-d356-430a-b906-55abbda8a940
# ╠═52deaf05-3549-49b8-bfe6-0f7c18df44eb
# ╠═7aa17393-80f2-45e3-ab3d-7dfd4adc2ed0
# ╠═534bbe60-cc91-45f4-9e90-7b0d2e576d59
# ╠═5859126d-a870-49c1-b74f-7eb50e29f78b
# ╠═77869159-c023-4097-8f2f-7af7739f6796
# ╠═f8bc4655-f11e-4af9-ae45-d8fdabd99cd2
# ╠═3f48ca1f-8d73-41b4-8721-d19ba881e670
# ╠═85ba6341-92d3-43dd-a9a1-a1e72a60b0d4
# ╠═5f0f00c9-05cc-4185-aad9-4bb0eab2c1ce
# ╠═e935064e-5e18-46ae-919d-83c8ac733f78
# ╠═915af65f-287a-4a11-9bb6-2755bc571a5a
# ╠═807f31dc-3e3f-4d2c-ab4a-4eb150f9599b
# ╠═8bb14d5a-1096-4fa8-a505-23b031bb0611
# ╠═a9dd3f32-92df-4ddf-9f35-6a8f8b73b5be
# ╠═e922cb53-d0de-44b9-bb79-c20dfac4263d
# ╠═48d03eee-4972-443f-a0c9-9fdaafd73c73
# ╠═21ea6b12-e2cb-46df-aed9-0ec6cd2e0cf8
