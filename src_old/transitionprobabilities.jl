# Second order system
# TODO: inverse function for f2!!!!
function get_ξ(e::EulerMaruyama,sde::SDE_Oscillator1D,t,t0,x,v) # x - f₁(ξ) = 0
    x - v*(t-t0)
end

function _tp(sde::SDE_Oscillator1D,par,x,t,x0,t0; method = EulerMaruyama(), kwargs...) # transition probability for scalar problem
    f2, g2 = method(sde,par,x...,t,x0...,t0)
    σ² = g2^2*(t-t0)
    μ = x0[2] + f2*(t-t0)
    normal1D(μ,σ²,x[2])
end

function _TP(sde::SDE_Oscillator1D,par,x,t,x0,t0;  method = EulerMaruyama(), kwargs...)
    ξ = get_ξ(method,sde,t,t0,x...)
    _tp(sde,par,x,t,(ξ,x0[2]),t0; method = method, kwargs...)
end

