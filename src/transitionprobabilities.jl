## Base distributions
# Normal Distribution
function normal1D(x) # 1D standard normal distribuiton
    exp(-0.5*x^2)/sqrt(2π)
end

function normal1D(μ,σ,x) # 1D normal distribuiton
    normal1D((x-μ)/σ) / σ
end

function normal1D_σ2(μ,σ²,x) # 1D standard normal distribuiton
    exp(-0.5*((x-μ)^2)/σ²)/sqrt(2π*σ²)
end
# Transformation for a Milstein step: y = x + Rx²; where x ∼ N(μₓ, σₓ)
function get_milstein_xvals(R,y)
    D2 = 1 + 4*R*y # D²
    if D2 > zero(D2)
        D = sqrt(D2)
        _2b = one(D)/(2R)
        D_2b = D*_2b

        x1,x2 = (- D_2b - _2b, D_2b - _2b)

        d = one(D)/(D)
        return x1, x2, d
    end
    
    z = zero(D2)
    return z, z, z
end

function milstein_dist(μₓ, σₓ, R, y)#, atol = sqrt(eps(R)))
    if abs2(R) < eps(R)
        return normal1D(μₓ, σₓ, y)
    end
    
    x1, x2, d = get_milstein_xvals(R,y)
    (normal1D(μₓ, σₓ, x1)  + normal1D(μₓ, σₓ, x2)) * d
end

## Step distributions
function _tp(sde::AbstractSDE{N,k},x,t,x0,t0; kwargs...) where {N,k} # transition probability for scalar problem
    _tp(sde,_par(sde),x,t,x0,t0; kwargs...)
end

## Scalar problems
function _tp(sde::SDE{1,1},par,x₁,t₁,x₀,t₀; method = EulerMaruyama()) # transition probability for scalar problem
    μ₀, σ₀ = method(sde,par,x₀,t₀,t₁-t₀)
    σ² = σ₀^2*(t₁-t₀)
    μ = x₀ + μ₀*(t₁-t₀)
    normal1D_σ2(μ,σ²,x₁)
end

# Second order system

function _tp(sde::SDE_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀; method = EulerMaruyama(), kwargs...)
    f₂, g₂ = method(sde,par,x₁,v₁,t₁,x₀,v₀,t₀)
    σ² = g₂^2*(t₁-t₀)
    μ = v₀ + f₂*(t₁-t₀)
    normal1D_σ2(μ,σ²,v₁)
end

# transition probability without impact/at wall
function _tp(sde::SDE_VI_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀, _r; method = EulerMaruyama(), Q_atwall::Bool = false, kwargs...)
    f₂, g₂ = method(sde.osc1D,par,x₁,v₁,t₁,x₀,v₀,t₀)
    Δt = (t₁-t₀);
    σ² = g₂^2*Δt
    # μ = v₁ + f₂*(t₁-t₀)
    μᵥ = v₀ + f₂*Δt
    res = normal1D_σ2(μᵥ,σ²,v₁)
    if Q_atwall
        r = _r(v₀) 
        res = res + normal1D(-r*μᵥ,σ²*r^2,v₁)/r
    end
    return res
end

# transition probability at impact
function _tp(sde::SDE_VI_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀, Δt₁, Δt₂, _r; method = EulerMaruyama(), kwargs...)
    r = _r(v₀)
    f₂, g₂ = method(sde,par,x₁,v₁,t₁,x₀,v₀,t₀,Δt₁,Δt₂,r)
    σ² = g₂^2*(r^2*Δt₁+Δt₂)
    μᵥ = -r*v₀ + f₂*(Δt₂-r*Δt₁)
    return normal1D_σ2(μᵥ,σ²,v₁) # abs(r) r should be positive!
end
