"""
 Estimators using local anisotropies
 - Kriging using non-stationary covariances: MW, KC, DE(?)
 - IDW: MW, averaged matrix
"""

## Necessary custom functions in the end
## Think how to deal with local pars for vaio structure


@estimsolver LocalKriging begin
  @param variogram = GaussianVariogram()
  @param mean = nothing #0.0 or [1.0,0.8,....]
  @param method = :MovingWindows
  @param localpars = [nothing]
  @param minneighbors = 1
  @param maxneighbors = 40
  @param neighborhood = nothing
  @param distance = Euclidean()
end

function solve(problem::EstimationProblem, solver::LocalKriging)
  # preprocess user input
  preproc = local_preprocess(problem, solver)

  # results for each variable
  μs = []; σs = []
  for (var, V) in variables(problem)
    varμ, varσ = local_solve_approx(problem, var, preproc)
    push!(μs, var => varμ)
    push!(σs, var => varσ)
  end

  EstimationSolution(domain(problem), Dict(μs), Dict(σs))
end


function local_preprocess(problem::EstimationProblem, solver::LocalKriging)
  # retrieve problem info
  pdomain = domain(problem)
  pdata = data(problem)

  # result of preprocessing
  preproc = Dict{Symbol,NamedTuple}()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine which Kriging variant to use
      if varparams.mean ≠ nothing
        estimator = SimpleKriging(varparams.variogram, varparams.mean)
      else
        estimator = OrdinaryKriging(varparams.variogram)
      end

      # determine minimum/maximum number of neighbors
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors

      # determine neighborhood search method
      if varparams.maxneighbors ≠ nothing
        if varparams.neighborhood ≠ nothing
          # local search with a neighborhood
          neigh = varparams.neighborhood

          if neigh isa BallNeighborhood
            bsearcher = KBallSearcher(pdata, maxneighbors, neigh)
          else
            searcher  = NeighborhoodSearcher(pdata, neigh)
            bsearcher = BoundedSearcher(searcher, maxneighbors)
          end
        else
          # nearest neighbor search with a distance
          distance = varparams.distance
          bsearcher = KNearestSearcher(pdata, maxneighbors, metric=distance)
        end
      else
        # use all data points as neighbors
        bsearcher = nothing
      end

      # local inputs
      method = varparams.method
      localpars = varparams.localpars

      # check pars
      okmeth = method in [:MovingWindows, :KernelConvolution]
      @assert okmeth "method must be :MovingWindows or :KernelConvolution"

      oklocal1 = length(localpars) == nelms(pdomain)
      oklocal2 = typeof(localpars[1]) <: LocalParameters
      @assert oklocal1 "number of local parameters must match domain points"
      @assert oklocal2 "wrong format of local parameters"

      # save preprocessed input
      preproc[var] = (estimator=estimator,
                      minneighbors=minneighbors,
                      maxneighbors=maxneighbors,
                      bsearcher=bsearcher,
                      method=method,
                      localpars=localpars)
    end
  end

  preproc
end


function local_solve_approx(problem::EstimationProblem, var::Symbol, preproc)
    # retrieve problem info
    pdata = data(problem)
    pdomain = domain(problem)
    N = ncoords(pdomain)
    T = coordtype(pdomain)

    # unpack preprocessed parameters
    estimator, minneighbors, maxneighbors, bsearcher, method, localpars = preproc[var]
    KC = method == :KernelConvolution ? true : false

    # if KC, pass localpars to hard data
    KC && (hdlocalpars = nnlocalpars(pdata,pdomain,localpars))

    # determine value type
    V = variables(problem)[var]

    # pre-allocate memory for result
    varμ = Vector{V}(undef, nelms(pdomain))
    varσ = Vector{V}(undef, nelms(pdomain))

    # pre-allocate memory for coordinates
    xₒ = MVector{N,T}(undef)

    # pre-allocate memory for neighbors
    neighbors = Vector{Int}(undef, maxneighbors)
    X = Matrix{T}(undef, N, maxneighbors)

    # estimation loop
    for location in traverse(pdomain, LinearPath())
      # coordinates of neighborhood center
      coordinates!(xₒ, pdomain, location)

      # find neighbors with previously estimated values
      nneigh = search!(neighbors, xₒ, bsearcher)

      localestimator = KC ? nothing : mwvario(estimator, localpars[location])

      # skip location in there are too few neighbors
      if nneigh < minneighbors
        varμ[location] = NaN
        varσ[location] = NaN
      else
        # final set of neighbors
        nview = view(neighbors, 1:nneigh)

        # get neighbors coordinates and values
        coordinates!(X, pdata, nview)

        Xview = view(X,:,1:nneigh)
        zview = view(pdata[var], nview)

        # fit estimator to data and predict mean and variance
        if !KC
          krig = local_fit(localestimator, Xview, zview)
          μ, σ² = local_predict(krig, xₒ)
        else
          println("To do")
          #localargs = view(hdlocalpars, nneighids)
          #push!(localargs, localpars[location])
          #krig = local_fit(estimator, Xview, zview, localargs)
          #μ, σ² = local_predict(krig, xₒ, localargs)
        end

        varμ[location] = μ
        varσ[location] = σ²
      end
    end

    varμ, varσ
end


function local_fit(estimator::KrigingEstimator, X::AbstractMatrix,
  z::AbstractVector, localpars::AbstractVector=[nothing])

  # build Kriging system
  LHS = local_lhs(estimator, X, localpars)
  RHS = Vector{eltype(LHS)}(undef, size(LHS,1))

  # factorize LHS
  FLHS = factorize(estimator, LHS)

  # record Kriging state
  state = KrigingState(X, z, FLHS, RHS)

  # return fitted estimator
  FittedKriging(estimator, state)
end


function local_lhs(estimator::KrigingEstimator, X::AbstractMatrix,
  localpars::AbstractVector)

  γ = estimator.γ
  nobs = size(X, 2)
  ncons = nconstraints(estimator)

  # pre-allocate memory for LHS
  x = view(X,:,1)
  T = Variography.result_type(γ, x, x)
  m = nobs + ncons
  LHS = Matrix{T}(undef, m, m)

  # set variogram/covariance block
  localpars[1] == nothing && Variography.pairwise!(LHS, γ, X)
  localpars[1] != nothing && kcfill!(LHS, γ, X, localpars)

  if isstationary(γ)
    for j=1:nobs, i=1:nobs
      @inbounds LHS[i,j] = sill(γ) - LHS[i,j]
    end
  end

  # set blocks of constraints
  set_constraints_lhs!(estimator, LHS, X)

  LHS
end


function set_local_rhs!(estimator::FittedKriging, xₒ::AbstractVector,
  localpars::AbstractVector)

  γ = estimator.estimator.γ
  X = estimator.state.X
  RHS = estimator.state.RHS

  #############################################################
  #### FILL SOME OTHER WAY FOR KC
  #############################################################

  # RHS variogram/covariance
  @inbounds for j in 1:size(X, 2)
    xj = view(X,:,j)
    RHS[j] = isstationary(γ) ? sill(γ) - γ(xj, xₒ) : γ(xj, xₒ)
  end

  set_constraints_rhs!(estimator, xₒ)
end

function local_weights(estimator::FittedKriging, xₒ::AbstractVector,
  localpars::AbstractVector)

  nobs = size(estimator.state.X, 2)

  set_local_rhs!(estimator, xₒ, localpars)

  # solve Kriging system
  x = estimator.state.LHS \ estimator.state.RHS

  λ = view(x,1:nobs)
  ν = view(x,nobs+1:length(x))

  KrigingWeights(λ, ν)
end

local_predict(estimator::FittedKriging, xₒ::AbstractVector, localpars::AbstractVector=[nothing]) =
  combine(estimator, local_weights(estimator, xₒ, localpars), estimator.state.z)



function nnlocalpars(pdata,pdomain,localpars)
  nothing
end
