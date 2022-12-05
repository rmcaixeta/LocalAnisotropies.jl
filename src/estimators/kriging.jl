# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from KrigingEstimators.jl
# ------------------------------------------------------------------

function solve(problem::EstimationProblem, solver::LocalKriging)
  # preprocess user input
  preproc = local_preprocess(problem, solver)

  # results for each variable
  μs = []; σs = []
  for var in name.(variables(problem))
    varμ, varσ = local_solve_approx(problem, var, preproc)
    push!(μs, var => varμ)
    push!(σs, Symbol(var,"-variance") => varσ)
  end

  georef((; μs..., σs...), domain(problem))
end


function local_preprocess(problem::EstimationProblem, solver::LocalKriging)
  # retrieve problem info
  pdomain = domain(problem)
  pdata   = data(problem)
  ndata   = nvals(pdata)

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
        # upper bound in maxneighbors
        maxneighbors > ndata && (maxneighbors = ndata)
        if varparams.neighborhood ≠ nothing
          # local search with a neighborhood
          neigh = varparams.neighborhood

          if neigh isa MetricBall
            bsearcher = KBallSearch(pdata, maxneighbors, neigh)
          else
            searcher  = BallSearch(pdata, neigh)
            bsearcher = BoundedSearch(searcher, maxneighbors)
          end
        else
          # nearest neighbor search with a distance
          distance = varparams.distance
          bsearcher = KNearestSearch(pdomain, maxneighbors, metric=distance)
        end
      else
        # use all data points as neighbors
        bsearcher = nothing
      end

      # local inputs
      method = varparams.method
      KC = method == :KernelConvolution ? true : false
      localaniso, localanisohd = (varparams.localaniso, varparams.localanisohd)
      (localanisohd != nothing && KC) && (localanisohd = toqmat(localanisohd))

      # check pars
      okmeth = method in [:MovingWindows, :KernelConvolution]
      @assert okmeth "method must be :MovingWindows or :KernelConvolution"

      oklocal1 = length(localaniso.rotation) == nvals(pdomain)
      oklocal2 = typeof(localaniso) <: LocalAnisotropy
      @assert oklocal1 "number of local anisotropies must match domain points"
      @assert oklocal2 "wrong format of local anisotropies"

      # save preprocessed input
      preproc[var] = (estimator=estimator,
                      minneighbors=minneighbors,
                      maxneighbors=maxneighbors,
                      bsearcher=bsearcher,
                      method=method,
                      localaniso=localaniso,
                      localanisohd=localanisohd)
    end
  end

  preproc
end


function local_solve_approx(problem::EstimationProblem, var::Symbol, preproc)
    # retrieve problem info
    pdata = data(problem)
    pdomain = domain(problem)

    mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

    # unpack preprocessed parameters
    estimator, minneighbors, maxneighbors, bsearcher, method, localaniso, hdlocalaniso = preproc[var]
    KC = method == :KernelConvolution ? true : false

    # if KC, pass localaniso to hard data
    (KC && hdlocalaniso==nothing) && (hdlocalaniso = grid2hd_qmat(pdata,pdomain,localaniso))

    # determine value type
    V = mactypeof[var]

    # pre-allocate memory for result
    varμ = Vector{V}(undef, nvals(pdomain))
    varσ = Vector{V}(undef, nvals(pdomain))

    # estimation loop
    Threads.@threads for location in traverse(pdomain, LinearPath())
      # pre-allocate memory
      neighbors = Vector{Int}(undef, maxneighbors)

      # centroid of neighborhood center
      pₒ = centroid(pdomain, location)

      # find neighbors with previously estimated values
      nneigh = search!(neighbors, pₒ, bsearcher)

      localpar = (rotation(localaniso,location),magnitude(localaniso,location))
      localestimator = KC ? nothing : mwvario(estimator, localpar)

      # skip location in there are too few neighbors
      if nneigh < minneighbors
        varμ[location] = missing
        varσ[location] = missing
      else
        # final set of neighbors
        nview = view(neighbors, 1:nneigh)

        # get neighbors centroid and values
        X = view(pdata, nview)
        println("NVALS: $(nvals(pdata)) $nneigh $neighbors")

        # not using block kriging yet, need more tests
        # uₒ = pdomain[location]

        # fit estimator to data and predict mean and variance
        if !KC
          krig = local_fit(localestimator, X)
          μ, σ² = local_predict(krig, var, pₒ)
        else
          ∑neighs = view(hdlocalaniso, nview)
          krig = local_fit(estimator, X, ∑neighs)
          Qx₀ = qmat(localpar...)
          μ, σ² = local_predict(krig, var, pₒ, (Qx₀,∑neighs))
        end

        varμ[location] = μ
        varσ[location] = σ²
      end
    end

    varμ, varσ
end


function local_lhs(estimator::KrigingEstimator, domain,
  localaniso::AbstractVector)

  γ = estimator.γ
  nobs = nvals(domain)
  ncons = nconstraints(estimator)

  # pre-allocate memory for LHS
  u = first(domain)
  T = Variography.result_type(γ, u, u)
  m = nobs + ncons
  println("Check: $nobs $ncons $m $(length(domain))")
  LHS = Matrix{T}(undef, m, m)

  # set variogram/covariance block
  KC = (localaniso[1] == nothing) ? false : true

  if !KC
    Variography.pairwise!(LHS, γ, domain)
    if isstationary(γ)
      for j=1:nobs, i=1:nobs
        @inbounds LHS[i,j] = sill(γ) - LHS[i,j]
      end
    end
  else
    kcfill!(LHS, γ, domain, localaniso)
  end

  # set blocks of constraints
  set_constraints_lhs!(estimator, LHS, domain)

  LHS
end


function set_local_rhs!(estimator::FittedKriging, pₒ,
  localaniso::Tuple)

  γ = estimator.estimator.γ
  X = domain(estimator.state.data)
  RHS = estimator.state.RHS

  KC = localaniso[1]==nothing ? false : true

  # RHS variogram/covariance
  if !KC
    @inbounds for j in 1:nvals(X)
      xj = centroid(X, j)
      RHS[j] = isstationary(γ) ? sill(γ) - γ(xj, pₒ) : γ(xj, pₒ)
    end
  else
    @inbounds for j in 1:nvals(X)
      xj = centroid(X, j)
      RHS[j] = kccov(γ, pₒ, xj, localaniso[1], localaniso[2][j])
    end
  end

  set_constraints_rhs!(estimator, pₒ)
end


function local_fit(estimator::KrigingEstimator, data,
  localaniso::AbstractVector=[nothing])

  D = domain(data)

  # build Kriging system
  LHS = local_lhs(estimator, D, localaniso)
  RHS = Vector{eltype(LHS)}(undef, size(LHS,1))

  # factorize LHS
  FLHS = factorize(estimator, LHS)

  # record Kriging state
  VARTYPE = Variography.result_type(estimator.γ, first(D), first(D))
  state = KrigingState(data, FLHS, RHS, VARTYPE)

  # return fitted estimator
  FittedKriging(estimator, state)
end


function local_weights(estimator::FittedKriging, pₒ,
  localaniso::Tuple)

  nobs = nvals(estimator.state.data)

  set_local_rhs!(estimator, pₒ, localaniso)

  # solve Kriging system
  x = estimator.state.LHS \ estimator.state.RHS

  λ = view(x,1:nobs)
  ν = view(x,nobs+1:length(x))

  KrigingWeights(λ, ν)
end


function local_predict(estimator::FittedKriging, var, pₒ, localaniso::Tuple=(nothing,))
  data = estimator.state.data
  combine(estimator, local_weights(estimator, pₒ, localaniso), data[var])
end
