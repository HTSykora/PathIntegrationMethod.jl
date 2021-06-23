# time-step methods
abstract type SDEMethod end
struct EulerMaruyama <: SDEMethod
end


# scalar
function normal1D(μ,σ²,x)
    exp(-0.5*((x - μ)^2 / σ²))/sqrt(2π*σ²)
end

function (e::EulerMaruyama)(sde::SDE{1,1},par,x0,t0)
    sde.f(x0,par,t0), sde.g(x0,par,t0)^2
end

function _tp(sde::SDE{1,1},par,x,t,x0,t0; method = EulerMaruyama()) # transition probability for scalar problem
    μ₀, σ₀² = method(sde,par,x0,t0)
    σ² = σ₀²*(t-t0)
    μ = x0 + μ₀*(t-t0)
    normal1D(μ,σ²,x)
end
function _tp(sde::SDE{1,1},x,t,x0,t0; kwargs...) # transition probability for scalar problem
    _tp(sde,sde.par,x,t,x0,t0; kwargs...)
end


# Second order system
# TODO: inverse function for f2!!!!
function (e::EulerMaruyama)(sde::SDE{2,2},par,x0,t0,t1)
    f1dt = x0[2]*(t1-t0)
    _x0 = copy(x0); _x0[1] = x0[1] - f1dt
    f2 = sde.f(2,_x0,par,t0)
    g2 = sde.g(2,_x0,par,t0)^2
    f2, g2
end

function _tp(sde::SDE{2,2},par,x,t,x0,t0; method = EulerMaruyama()) # transition probability for scalar problem
    f2, g2 = method(sde,par,x0,t0,t1)
    σ² = g2*(t-t0)
    μ = x0[2] + f2
    normal1D(μ,σ²,x[2])
end



