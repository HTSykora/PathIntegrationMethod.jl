function normal1D(μ,σ²,x) # 1D Gaussian distribuiton
    exp(-0.5*((x - μ)^2 / σ²))/sqrt(2π*σ²)
end

function _tp(sde::AbstractSDE{N,k},x,t,x0,t0; kwargs...) where {N,k} # transition probability for scalar problem
    _tp(sde,sde.par,x,t,x0,t0; kwargs...)
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

function _tp(sde::SDE_VI_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀, Q_impact, Δt₁, Δt₂, r; method = EulerMaruyama(), kwargs...)
    if Q_impact
        f₂, g₂ = method(sde,par,x₁,v₁,t₁,x₀,v₀,t₀,Δt₁,Δt₂,r)
        σ² = g₂^2*(r^2*Δt₁+Δt₂)
        μ = v₀ + f₂*(r*Δt₁+Δt₂)
        return normal1D(μ,σ²,v₁)/r # abs(r) r should be positive!
    else
        f₂, g₂ = method(sde.osc1D,par,x₁,v₁,t₁,x₀,v₀,t₀)
        σ² = g₂^2*(t₁-t₀)
        μ = v₀ + f₂*(t₁-t₀)
        return normal1D(μ,σ²,v₁)
    end
end
