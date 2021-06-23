# scalar

function _tp(sde::SDE{1,1},par,x,t,x0,t0; method = EulerMaruyama()) # transition probability for scalar problem
    μ₀, σ₀² = method(sde,par,x0,t0)
    σ² = σ₀²*(t-t0)
    μ = x0 + μ₀*(t-t0)
    exp(-0.5*((x - μ)^2 / σ²))/sqrt(2π*σ²)
end
function _tp(sde::SDE{1,1},x,t,x0,t0; kwargs...) # transition probability for scalar problem
    _tp(sde,sde.par,x,t,x0,t0; kwargs...)
end

# time-step methods
abstract type SDEMethod end
struct EulerMaruyama <: SDEMethod
end

function (e::EulerMaruyama)(sde::SDE{1,1},par,x0,t0)
    sde.f(x0,par,t0), sde.g(x0,par,t0)^2
end