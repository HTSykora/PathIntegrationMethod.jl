function normal1D(μ,σ²,x) # 1D Gaussian distribuiton
    exp(-0.5*((x - μ)^2 / σ²))/sqrt(2π*σ²)
end

function _tp(sde::AbstractSDE{N,k},x,t,x0,t0; kwargs...) where {N,k} # transition probability for scalar problem
    _tp(sde,sde.par,x,t,x0,t0; kwargs...)
end
function _TP(sde::AbstractSDE{N,k},x,t,x0,t0; kwargs...) where {N,k}
    _TP(sde,sde.par,x,t,x0,t0; kwargs...)
end

## Scalar problems
function _tp(sde::SDE{1,1},par,x,t,x0,t0; method = EulerMaruyama()) # transition probability for scalar problem
    μ₀, σ₀ = method(sde,par,x0,t0,t-t0)
    σ² = σ₀^2*(t-t0)
    μ = x0 + μ₀*(t-t0)
    normal1D(μ,σ²,x)
end

# Second order system
# TODO: inverse function for f2!!!!
function get_ξ(e::EulerMaruyama,sde::SDE_Oscillator1D,t₁,t₀,x₁,v₁,x₀,v₀) # x - f₁(ξ) = 0
    x₁ - v₀*(t₁-t₀)
end

function _tp(sde::SDE_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀; method = EulerMaruyama(), kwargs...) # transition probability for scalar problem
    f₂, g₂ = method(sde,par,x₁,v₁,t₁,x₀,v₀,t₀)
    σ² = g₂^2*(t₁-t₀)
    # μ = v₁ + f₂*(t₁-t₀)
    μ = v₀ + f₂*(t₁-t₀)
    normal1D(μ,σ²,v₁)
end

# function _TP(sde::SDE_Oscillator1D,par,x,t,x0,t0;  method = EulerMaruyama(), kwargs...)
#     ξ = get_ξ(method,sde,t,t0,x...)
#     _tp(sde,par,x,t,(ξ,x0[2]),t0; method = method, kwargs...)
# end
