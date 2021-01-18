using LocalAnisotropies
using GeoStatsBase
using Test

@testset "LocalAnisotropies.jl" begin
    ## convert data to LocalParameters
    dummy = georef((az=1:10, r1=1:10, r2=1:10), PointSet(rand(2,10)))
    pars  = LocalParameters(dummy, [:az], [:r1,:r2], :GSLIB)
    angs  = convertangles([30,30,30], :GSLIB, :Datamine)
    angs_ = convertangles.(pars.rotation, :GSLIB)
    # exportpars("C:\\Users\\test.csv", pars, :GSLIB)

    ## 2-D data
    # georef rand 2d data
    # sample 2d data
    # testset

    ## 3-D data
    # georef rand 3d data
    # sample 3d data
    # empty downscaled grid
    # localparameters(Grad)
    # localparameters(Geometric)
    # rescale_magnitude()
    # LocalKriging(MW)
    # LocalKriging(KC)
    # deformspace(M1)
    # deformspace(M2)
    # deformspace(M3)
    # Kriging(M1)
    # Kriging(M2)
    # Kriging(M3)

    ## interpolate pars
    # idwpars()
    # smoothpars()

end
