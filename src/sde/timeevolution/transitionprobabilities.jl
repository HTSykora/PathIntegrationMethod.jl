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

# Transitionprobability of a step
# function transitionprobability(step::SDEStep{d,k,m,sdeT,method}, x::Vararg{Any,N}) where{d, m, sdeT, method<:DiscreteTimeStepping{TDrift, TDiff}, N} where {TDrift, TDiff<:Maruyama}
#     # TODO: case with diagonal noise
#     # N == d-k+1
#     # σ2 = ...
#     # prod(normal1D_σ2(step.x1[k+i-1],σ2[i,i],x[i]) for i in 1:N)

#     # TODO: case with general noise
# end

# Single noise source only on the last coordinate
function transitionprobability(step::SDEStep{d,d,m,sdeT,method},x) where {d,m,sdeT, method<:DiscreteTimeStepping{TDrift, TDiff}} where {TDrift, TDiff<:Maruyama}
    σ2 = _Δt(step) * (get_g(step.sde)(d, step.x0,_par(step),_t0(step))^ 2)
    # detJ_correction = _detJ(step, x)
    normal1D_σ2(step.x1[d], σ2, x[d]);
end