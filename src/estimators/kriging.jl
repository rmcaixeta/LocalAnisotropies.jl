# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsModels.jl
# ------------------------------------------------------------------


struct LocalKrigingModel <: GeoStatsModel
    method::Symbol
    localaniso::LocalAnisotropy
    γ::Variogram
    skmean::Union{Number,Nothing}  # for later: local mean array
    hdlocalaniso::Union{AbstractVector,Nothing}
end

"""
    LocalKriging(params ...)

LocalKriging estimation solver, where `var::Symbol` is the variable name
and `param` is a `NamedTuple` containing the parameters below:

## Parameters

* `variogram` - Reference variogram model
* `mean`      - Simple Kriging mean
* `method`    - LocalKriging method. :MovingWindows or :KernelConvolution
  (default to :MovingWindows)
* `localaniso`    - Local parameters of the domain
* `localanisohd`  - Local parameters of the samples. Only necessary for
  :KernelConvolution method. They are automatically passed via NN from
  `localaniso` if not informed.
"""
LocalKriging(
    method::Symbol,
    localaniso::LocalAnisotropy,
    γ::Variogram;
    skmean = nothing,
    hdlocalaniso = nothing,
) = LocalKrigingModel(method, localaniso, γ, skmean, hdlocalaniso)


function nconstraints(::LocalKrigingModel)
    n = skmean == nothing ? 1 : 0 # OK otherwise SK
    n
end

function local_fit(model_::LocalKrigingModel, data; i, m)
    MW = (model_.method == :MovingWindows)
    localaniso = model_.localaniso
    localpar = (rotation(localaniso, i), magnitude(localaniso, i))
    if MW
        model = mwvario(model_, localpar)
    else
        model =
            model_.skmean == nothing ? GeoStatsModels.OrdinaryKriging(model_.γ) :
            GeoStatsModels.SimpleKriging(model_.γ, model_.skmean)
        hdlocalaniso = view(model_.hdlocalaniso, m)
        Qx₀ = qmat(localpar...)
    end

    γ = model.γ
    D = domain(data)

    # build Kriging system
    LHS = MW ? lhs(model, D) : local_lhs(model, D, hdlocalaniso)
    RHS = Vector{eltype(LHS)}(undef, size(LHS, 1))

    # factorize LHS
    FLHS = GeoStatsModels.factorize(model, LHS)

    # variance type
    VARTYPE = GeoStatsFunctions.returntype(γ, first(D), first(D))

    # record Kriging state
    state = KrigingState(data, FLHS, RHS, VARTYPE)

    # return fitted model
    MW ? FittedKriging(model, state) : LocalFittedKriging(model, state, Qx₀, hdlocalaniso)
end

struct LocalFittedKriging#{M<:LocalKrigingModel,S<:KrigingState}
    model::KrigingModel
    state::KrigingState
    Qx₀::Any
    hdlocalaniso::Any
end

FKC(m::LocalFittedKriging) = FittedKriging(m.model, m.state)

status(fitted::LocalFittedKriging) = issuccess(fitted.state.LHS)

predict(fitted::LocalFittedKriging, var, uₒ) =
    predictmean(FKC(fitted), weights(fitted, uₒ), var)

function predictprob(fitted::LocalFittedKriging, var, uₒ)
    w = local_weights(fitted, uₒ)
    μ = predictmean(fitted, w, var)
    σ² = predictvar(fitted, w)
    Normal(μ, √σ²)
end

predictvar(fitted::LocalFittedKriging, weights::KrigingWeights) =
    GeoStatsModels.predictvar(FKC(fitted), weights::KrigingWeights)

predictmean(fitted::LocalFittedKriging, weights::KrigingWeights, var) =
    GeoStatsModels.predictmean(FKC(fitted), weights::KrigingWeights, var)

set_constraints_rhs!(fitted::LocalFittedKriging, pₒ) =
    GeoStatsModels.set_constraints_rhs!(FKC(fitted), pₒ)

function local_lhs(model::KrigingModel, domain, localaniso)
    γ = model.γ
    nobs = nvals(domain)
    ncons = nconstraints(model)

    # pre-allocate memory for LHS
    u = first(domain)
    T = GeoStatsFunctions.returntype(γ, u, u)
    m = nobs + ncons
    LHS = Matrix{T}(undef, m, m)

    # set variogram/covariance block
    kcfill!(LHS, γ, domain, localaniso)

    # set blocks of constraints
    GeoStatsModels.set_constraints_lhs!(model, LHS, domain)

    LHS
end

function set_local_rhs!(fitted::LocalFittedKriging, pₒ)
    localaniso = (fitted.Qx₀, fitted.hdlocalaniso)
    γ = fitted.model.γ
    X = domain(fitted.state.data)
    RHS = fitted.state.RHS

    # RHS variogram/covariance
    @inbounds for j = 1:nvals(X)
        xj = centroid(X, j)
        RHS[j] = kccov(γ, pₒ, xj, localaniso[1], localaniso[2][j])
    end

    set_constraints_rhs!(fitted, pₒ)
end


function weights(fitted::LocalFittedKriging, pₒ)
    localaniso = (fitted.Qx₀, fitted.hdlocalaniso)
    nobs = nvals(fitted.state.data)

    set_local_rhs!(fitted, pₒ)

    # solve Kriging system
    x = fitted.state.LHS \ fitted.state.RHS

    λ = view(x, 1:nobs)
    ν = view(x, nobs+1:length(x))

    KrigingWeights(λ, ν)
end
