"""
 Extract local anisotropy from some reference input
 - Gradients 2D and 3D
 - PCA for scattered points

 Store parameters as quaternions + ratios

- Need to resolve:
- best design to scale ratios easily?
    fix reference vario and scale between ]0,1]?
    fix vario range to one and scale between [min,max] ranges?
- how to assign local pars differently to each structures? mayber better to only
allow if the ratios should rescale all or given list of structures


"""

struct LocalParameters
  rotation::AbstractVector{Quaternion}
  magnitude::AbstractArray{AbstractFloat,2}
end

slicelp(lp::LocalParameters,ids) = LocalParameters(lp.rotation[ids],lp.magnitude[:,ids])


function gridneighbors(img, i::CartesianIndex, window::Int)
    dim = Size(img)
    minid(d) = max(1,i[d]-window)
    maxid(d) = min(dim[d],i[d]+window)
    idx = Tuple([minid(d):maxid(d) for d in 1:length(dim)])
    CartesianIndices(idx)
end

function gradients(preimg, prop, window)
    # get dimensions
    dims = preimg.domain.dims
    N = length(dims)

    # extract gradients (maybe need an extra method to deal with big datasets)
    img = reshape(preimg.table[prop], Size(dims))
    g = imgradients(img, KernelFactors.sobel, "replicate")

    quat = Array{Quaternion}(undef,size(img)) # make some better way to store it
    m = Array{Vector}(undef,size(img)) # make some better way to store it

    Threads.@threads for i in CartesianIndices(img)
        tensor = zeros(Float64,N,N)
        for ng in gridneighbors(img, i, window)
            for j in CartesianIndices((1:N,1:N))
                tensor[j] += g[j[1]][ng]*g[j[2]][ng]
            end
        end

        T = eigen(SMatrix{N,N}(tensor))
        V = T.vectors[:, sortperm(T.values)]'
        if N==3
            eigv = V
        else
            eigv = zeros(Float64,3,3)
            eigv[1:2,1:2] = V
            eigv[3,3] = 1.0
            eigv = SMatrix{3,3}(eigv)
        end
        q = dcm_to_quat(eigv)
        λ = sort(T.values, rev=true)
        quat[i] = q
        m[i] = λ/λ[1]
    end

    LocalParameters(vec(quat), reduce(hcat,vec(m)))
end


# function geometry(ref_pt, pars)
#     # PCA on mesh points
# end


function rescale_magnitude(lp::LocalParameters, bounds1, bounds2=nothing)
    N = bounds2==nothing ? 2 : size(lp.magnitude,1)
    m = lp.magnitude
    bounds = [bounds1,bounds2]

    for i in 2:N
        mi = m[i,:]
        qx = quantile(mi, [0.05,0.95])
        mi[mi .< qx[1]] .= qx[1]
        mi[mi .> qx[2]] .= qx[2]
        r = (mi .- qx[1]) ./ (qx[2]-qx[1])
        typeof(bounds[i-1]) <: Number && (bounds[i-1]=(bounds[i-1],bounds[i-1]))
        b1, b2 = bounds[i-1]
        m[i,:] .= (b1 .+ (b2-b1) .* r)
    end

    LocalParameters(lp.rotation,m)
end


#Generic function that cross-validate multiple scenarios and return an optimal local field

Base.@kwdef mutable struct TestSet
  refimgs::AbstractVector # [(img,prop1),(img,prop2)]
  smooth::AbstractVector{Number}
  variogram::AbstractVector # [(:Y,γ1),(:X,γ2)]

  ratio1::AbstractArray = [(0.2,1.0),(0.5,1.0)]
  ratio2::AbstractArray = [nothing]
  meths::AbstractVector{Symbol} = [:MovingWindows,:KernelConvolution]
  maxneighbors::AbstractArray{Int} = [15,25,40]

  #search::AbstractArray{Search}=[KNN(15),KNN(30),KNN(50)]

  folds::Int = 10
  forcegradients::Bool = false

end


macro name(arg)
   string(arg)
end

function localpars(problem::EstimationProblem, t::TestSet)
    ip = Iterators.product
    vars = [v for (v,V) in variables(problem)]

    bestε = Dict(v => Inf for v in vars)
    optpars = Dict()
    listε = Dict(v => [] for v in vars)

    ids = grid2hd_ids(data(problem),domain(problem))

    for (img,w) in ip(t.refimgs, t.smooth)
        M = typeof(img[1]) <: RegularGrid || t.forcegradients
        #lpars = M ? gradients(img[1],img[2],w) : geometry(img,w)
        lpars = gradients(img[1],img[2],w)

        for (r1,r2) in ip(t.ratio1, t.ratio2)
            lx = rescale_magnitude(lpars, r1, r2)

            for (v,n,e,γ) in ip(vars, t.maxneighbors, t.meths, t.variogram)
                p = (variogram=γ, localpars=lx, method=e, maxneighbors=n)
                solver = LocalKriging(v => p)
                cv = CrossValidation(t.folds)
                ε = cverror(solver, problem, cv)
                push!(listε[v],ε[v])

                if bestε[v]>ε[v]
                    bestε[v] = ε[v]
                    optpars[v] = p
                end
            end
        end
    end

    optpars
end


## Interpolate LocalParameters: NN and IDW
