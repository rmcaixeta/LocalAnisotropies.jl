# ------------------------------------------------------------------
# Temporarily ported and adapted from GeoStatsBase.jl
# Licensed under the ISC License. See LICENSE in the project root.
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

  # error for a fold k
  function ε(k)
    # source and target indices
    sinds = [ind for i in [1:k-1; k+1:nfolds] for ind in folds[i]]
    tinds = folds[k]

    # source and target data
    train = view(sdata, sinds)
    hold  = view(sdata, tinds)

    # setup and solve sub-problem
    subproblem = EstimationProblem(train, domain(hold), Tuple(ovars))
    solution   = solve(subproblem, solver)

    # loss for each variable
    losses = map(ovars) do var
      y = hold[var]
      ŷ = EP ? solution[var].mean : solution[var]
      ℒ = value(loss[var], y, ŷ, AggMode.Mean())
      var => ℒ
    end

    Dict(losses)
  end

  # compute error for each fold in parallel
  εs = foldxt(vcat, Map(ε), 1:nfolds)

  # combine error from different folds
  Dict(var => mean(get.(εs, var, 0)) for var in ovars)
end
