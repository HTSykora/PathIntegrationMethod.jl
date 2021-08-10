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
function get_ξ(e::EulerMaruyama,sde::SDE_Oscillator1D,t₁,t₀,x₁,v₁,x₀,v₀) # x₁ - φₓ(ξ,v₀,t₀,t₁) = 0
    x₁ - v₀*(t₁-t₀)
end

function (e::EulerMaruyama)(sde::SDE_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀)
    f2 = sde.f([x₀,v₀],par,t₀)
    g2 = sde.g([x₀,v₀],par,t₀)
    f2, g2
end

## 1D VI Oscillator 
# x₁ - φₓ(ξ,v₀,t₀,t₁) = 0
function get_ξ(e::EulerMaruyama,sde::SDE_VI_Oscillator1D{wTT},t₁,t₀,x₁,v₁,x₀,v₀) where wTT# x - f₁(ξ) = 0
    ξ = x₁ - v₀*(t₁-t₀)
    Q_impact, wallID = hit(ξ,sde.w)
    if Q_impact
        d = sde.w[wallID].d;
        ξ = d - v₀*(t₁ - t₀) - (x₁ - d)
        Δt1 = sde.w[wallID].d - ξ
        Δt2 = 
    else
        Δt1 = nothing
        Δt2 = nothing
    end
    ξ, Q_impact, Δt1, Δt2
end
ξ, Q_impact, Δt1, Δt2 = get_ξ(IK.method,IK.sde,IK.t₁,IK.t₀,IK.pdgrid.xs[1][IK.idx₁[1]],nothing,nothing,v₀) # v₁, x₀ = nothing


function (e::EulerMaruyama)(sde::SDE_VI_Oscillator1D,par,x₁,v₁,t₁,x₀,v₀,t₀)
    f2 = sde.f([x₀,v₀],par,t₀)
    g2 = sde.g([x₀,v₀],par,t₀)
    f2, g2
end
