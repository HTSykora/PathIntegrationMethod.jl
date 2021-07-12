## Scalar problems
function (e::EulerMaruyama)(sde::SDE{1,1},par,x0,t0,Δt)
    sde.f(x0,par,t0), sde.g(x0,par,t0)
end
function (e::Milstein)(sde::SDE{1,1},par,x0,t0,Δt)
    f1 = sde.f(x0,par,t0)
    K = x0 + Δt*f1
    L = sde.g(x0,par,t0)
    sqΔt = sqrt(Δt)
    mil_corretion = sde.g(K + L*sqΔt,par,t0) - L
    f1, (L + mil_corretion*sqΔt)
end
function (e::RKMaruyama)(sde::SDE{1,1},par,x0,t0,Δt)
    k1 = sde.f(x0,par,t0)
    k2 = sde.f(x0 + 0.5*Δt*k1,par,t0+0.5*Δt)
    k3 = sde.f(x0 + 0.5*Δt*k2,par,t0+0.5*Δt)
    k4 = sde.f(x0 + Δt*k3,par,t0+Δt)
    
    (k1+k2+k3+k4)/6, sde.g(x0,par,t0)
end

## 1D Oscillator
function (e::EulerMaruyama)(sde::SDE_Oscillator1D,par,x,v,t,x0,v0,t0)
    f2 = sde.f([x0,v0],par,t0)
    g2 = sde.g([x0,v0],par,t0)
    f2, g2
end
