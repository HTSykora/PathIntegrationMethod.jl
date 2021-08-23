function normal1D(μ,σ²,x) # 1D Gaussian distribuiton
    exp(-0.5*((x - μ)^2 / σ²))/sqrt(2π*σ²)
end

function _tp(sde::AbstractSDE{N,k},x,t,x0,t0; kwargs...) where {N,k} # transition probability for scalar problem
    _tp(sde,_par(sde),x,t,x0,t0; kwargs...)
end

## Scalar problems
function _tp(sde::SDE{1,1},par,x₁,t₁,x₀,t₀; method = EulerMaruyama()) # transition probability for scalar problem
    μ₀, σ₀ = method(sde,par,x₀,t₀,t₁-t₀)
    σ² = σ₀^2*(t₁-t₀)
    μ = x₀ + μ₀*(t₁-t₀)
    normal1D(μ,σ²,x₁)
end

# Second order system

function _tp(sde::SDE_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀; method = EulerMaruyama(), kwargs...)
    f₂, g₂ = method(sde,par,x₁,v₁,t₁,x₀,v₀,t₀)
    σ² = g₂^2*(t₁-t₀)
    # μ = v₁ + f₂*(t₁-t₀)
    μ = v₀ + f₂*(t₁-t₀)
    normal1D(μ,σ²,v₁)
end

# transition propability without impact
function _tp(sde::SDE_VI_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀, _r; method = EulerMaruyama(), Q_atwall::Bool = false, kwargs...)
    f₂, g₂ = method(sde.osc1D,par,x₁,v₁,t₁,x₀,v₀,t₀)
    σ² = g₂^2*(t₁-t₀)
    # μ = v₁ + f₂*(t₁-t₀)
    μ = v₀ + f₂*(t₁-t₀)
    res = normal1D(μ,σ²,v₁)
    if Q_atwall
        r = _r(v₀) 
        return res + normal1D(-r*μ,σ²*r^2,v₁)
    else
        return res
    end
end
# transition propability at impact
function _tp(sde::SDE_VI_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀, Δt₁, Δt₂, _r; method = EulerMaruyama(), kwargs...)
    r = _r(v₀)
    f₂, g₂ = method(sde,par,x₁,v₁,t₁,x₀,v₀,t₀,Δt₁,Δt₂,r)
    σ² = g₂^2*(r^2*Δt₁+Δt₂)
    μ = -r*v₀ + f₂*(Δt₂-r*Δt₁)
    return normal1D(μ,σ²,v₁) # abs(r) r should be positive!
end
