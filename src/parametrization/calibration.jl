# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# Need to be refactored and improved

Base.@kwdef mutable struct TestPars
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

function localanisotropies(t::TestPars, problem::EstimationProblem)
    ip = Iterators.product
    vars = [v for (v,V) in variables(problem)]

    bestε = Dict(v => Inf for v in vars)
    optpars = Dict()
    listε = Dict(v => [] for v in vars)

    ids = grid2hd_ids(data(problem),domain(problem))

    for (img,w) in ip(t.refimgs, t.smooth)
        M = typeof(img[1]) <: CartesianGrid || t.forcegradients
        #lpars = M ? gradients(img[1],img[2],w) : geometry(img,w)
        lpars = localanisotropies(Gradients,img[1],img[2],w)

        for (r1,r2) in ip(t.ratio1, t.ratio2)
            lx = rescale_magnitude(lpars, r1, r2)

            for (v,n,e,γ) in ip(vars, t.maxneighbors, t.meths, t.variogram)
                p = (variogram=γ, localaniso=lx, method=e, maxneighbors=n)
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
