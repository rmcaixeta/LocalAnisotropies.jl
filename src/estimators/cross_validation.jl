# ------------------------------------------------------------------
# Temporarily ported and adapted from GeoStatsBase.jl
# ------------------------------------------------------------------


function cverror(solver::AbstractSolver,
               problem::EstimationProblem,
               eestimator::CrossValidation)
  # problem info
  sdata = data(problem)
  ovars = [v for (v,V) in variables(problem)]

  # retrieve problem info
  partitioner = eestimator.partitioner
  loss  = eestimator.loss
  for var in ovars
    if var ∉ keys(loss)
      loss[var] = defaultloss(sdata[var][1])
    end
  end

  # folds for cross-validation
  folds  = subsets(partition(sdata, partitioner))
  nfolds = length(folds)
  dataids = grid2hd_ids(sdata,domain(problem))
  hdpars = Dict(v => slice(solver.vparams[v].localpars, dataids) for v in ovars)

  # error for a fold k
  function ε(k)
    # source and target indices
    sinds = [ind for i in [1:k-1; k+1:nfolds] for ind in folds[i]]
    tinds = folds[k]

    # source and target data
    train = view(sdata, sinds)
    hold  = view(sdata, tinds)
    vpars = []
    for v in ovars
      p = pars2tuple(solver.vparams[v])
      # find one to modify solver; or create new one
      p = @set p.localparshd = slice(hdpars[v], sinds)
      p = @set p.localpars = slice(hdpars[v], tinds)
      push!(vpars,p)
    end

    lsolver = LocalKriging(zip(ovars,vpars)...)

    # setup and solve sub-problem
    subproblem = EstimationProblem(train, domain(hold), Tuple(ovars))
    solution   = solve(subproblem, lsolver)

    # loss for each variable
    losses = map(ovars) do var
      y = hold[var]
      ŷ = solution[var].mean
      ℒ = value(loss[var], y, ŷ, AggMode.Mean())
      var => ℒ
    end

    Dict(losses)
  end

  # compute error for each fold in parallel
  εs = mapreduce(ε, vcat, 1:nfolds)

  # combine error from different folds
  Dict(var => mean(get.(εs, var, 0)) for var in ovars)
end


pars2tuple(p) =
  (;(v=>getfield(p, v) for v in fieldnames(typeof(p)) if v != :__dummy__)...)
