"""
 Estimators using local anisotropies
 - Kriging using non-stationary covariances: MW, KC, DE(?)
 - IDW: MW, averaged matrix
"""


# Just change Mahalanobis locally?
      # Ok for MW
      # Need to add extra value for covar in KC



"""
# Loop locations
function solve_approx(problem::EstimationProblem, var::Symbol, preproc) end
function solve_exact(problem::EstimationProblem, var::Symbol, preproc) end

# Build matrices
function fit(estimator::KrigingEstimator, X::AbstractMatrix, z::AbstractVector)
      # build Kriging system
      LHS = lhs(estimator, X)
      RHS = Vector{eltype(LHS)}(undef, size(LHS,1))
end

function lhs(estimator::KrigingEstimator, X::AbstractMatrix)
  pairwise!(LHS, γ, X)
end
"""
