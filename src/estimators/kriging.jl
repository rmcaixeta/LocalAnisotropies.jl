# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# Adapted from GeoStatsModels.jl
# ------------------------------------------------------------------


abstract type LocalKrigingModel <: GeoStatsModel end

struct MW_OKModel{G<:Variogram} <: LocalKrigingModel
    localaniso::LocalAnisotropy
    γ::G
end

struct MW_SKModel{G<:Variogram,V} <: LocalKrigingModel
    localaniso::LocalAnisotropy
    γ::G
    μ::V
end

struct KC_OKModel{G<:Variogram} <: LocalKrigingModel
    localaniso::LocalAnisotropy
    γ::G
    hdlocalaniso::Union{AbstractVector,Nothing}
end

struct KC_SKModel{G<:Variogram,V} <: LocalKrigingModel
    localaniso::LocalAnisotropy
    γ::G
    μ::V
    hdlocalaniso::Union{AbstractVector,Nothing}
end

MWModels = Union{MW_OKModel,MW_SKModel}
KCModels = Union{KC_OKModel,KC_SKModel}

krig_estimator(model::MWModels, localpar) = mw_estimator(model, localpar)
krig_estimator(model::KC_OKModel, localpar = nothing) = OrdinaryKriging(model.γ)
krig_estimator(model::KC_SKModel, localpar = nothing) = SimpleKriging(model.γ, model.μ)

neighs_localaniso(model::MWModels, m) = nothing
neighs_localaniso(model::KCModels, m) = view(model.hdlocalaniso, m)

function neighs_localaniso(
    localaniso::LocalAnisotropy,
    dom,
    neigh;
    method::Symbol = :KernelConvolution,
)
    method != :KernelConvolution ? nothing : qmat(nnpars(localaniso, dom, neigh))
end

function initmodel(model::KCModels, geotable, pdomain)
    hd = grid2hd_qmat(geotable, pdomain, model.localaniso)
    la = model.localaniso
    γ = model.γ
    model isa KC_OKModel ? KC_OKModel(la, γ, hd) : MW_SKModel(la, γ, model.μ, hd)
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
function LocalKriging(
    method::Symbol,
    localaniso::LocalAnisotropy,
    γ::Variogram;
    μ = nothing,
    hdlocalaniso::Union{LocalAnisotropy,Nothing} = nothing,
)
    if method == :MovingWindows
        isnothing(μ) ? MW_OKModel(localaniso, γ) : MW_SKModel(localaniso, γ, μ)
    elseif method == :KernelConvolution
        hd = qmat(hdlocalaniso)
        isnothing(μ) ? KC_OKModel(localaniso, γ, hd) : MW_SKModel(localaniso, γ, μ, hd)
    else
        @assert false "method must be :MovingWindows or :KernelConvolution"
    end
end

nconstraints(::LocalKrigingModel) = isnothing(μ) ? 1 : 0 # OK otherwise SK


function local_fit(model_::LocalKrigingModel, data; i, m)
    localaniso = model_.localaniso
    localpar = localpair(localaniso, i)
    model = krig_estimator(model_, localpar)

    γ = model.γ
    D = domain(data)

    hdlocalaniso = neighs_localaniso(model_, m)
    Qx₀ = model_ isa KCModels ? qmat(localpar...) : nothing

    # build Kriging system
    LHS = local_lhs(model, D, hdlocalaniso)
    RHS = Vector{eltype(LHS)}(undef, size(LHS, 1))

    # factorize LHS
    FLHS = GeoStatsModels.factorize(model, LHS)

    # variance type
    VARTYPE = GeoStatsFunctions.returntype(γ, first(D), first(D))

    # record Kriging state
    state = KrigingState(data, FLHS, RHS, VARTYPE)

    # return fitted model
    localfitting(model, state, Qx₀, hdlocalaniso)
end

localfitting(model, state, q, hd) =
    isnothing(q) ? FittedKriging(model, state) : LocalFittedKriging(model, state, q, hd)

struct LocalFittedKriging{Q,H}
    model::KrigingModel
    state::KrigingState
    Qx₀::Q
    hdlocalaniso::H
end

FKC(m::LocalFittedKriging) = FittedKriging(m.model, m.state)

local_status(fitted::LocalFittedKriging) = issuccess(fitted.state.LHS)
local_status(fitted::FittedKriging) = GeoStatsModels.status(fitted)

predict(fitted::LocalFittedKriging, var, uₒ) =
    predictmean(FKC(fitted), weights(fitted, uₒ), var)

function predictprob(fitted::LocalFittedKriging, var, uₒ)
    w = weights(fitted, uₒ)
    μ = predictmean(fitted, w, var)
    σ² = predictvar(fitted, w)
    σ² = fitted.model isa OrdinaryKriging  ? σ² + 2 * w.ν[1] : σ²
    Normal(μ, √σ²)
end

predictvar(fitted::LocalFittedKriging, weights::KrigingWeights) =
    GeoStatsModels.predictvar(FKC(fitted), weights::KrigingWeights)

predictmean(fitted::LocalFittedKriging, weights::KrigingWeights, var) =
    GeoStatsModels.predictmean(FKC(fitted), weights::KrigingWeights, var)

set_constraints_rhs!(fitted::LocalFittedKriging, pₒ) =
    GeoStatsModels.set_constraints_rhs!(FKC(fitted), pₒ)


function local_lhs(model::KrigingModel, domain, localaniso)
    if isnothing(localaniso)
        lhs(model, domain)
    else
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
