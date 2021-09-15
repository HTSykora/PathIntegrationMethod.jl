function eval_driftstep!(step::SDEStep{d,k,m, sdeT, DiscreteTimeStepping{Euler}}) where {d,k,m,sdeT}
    Δt = _Δt(step)
    for i in 1:d
        step.x1[i] = step.sde.f(i,step.x0,_par(step),_t0(step))*Δt
    end
end


function eval_driftstep_xI_sym(sde::AbstractSDE{d,k,m}, method::DiscreteTimeStepping{<:Euler}, x, par, t0, t1) where {d,k,m}
    [x[i] + sde.f(i,x,par,t0)*(t1-t0) for i in 1:k-1]
end
function update_drift_x!(step::SDEStep{d,k,m, sdeT, DiscreteTimeStepping{Euler}}) where {d,k,m,sdeT}
    for i in 1:k-1
        step.x0[i] = step.steptracer.temp[i]
    end
    Δt = _Δt(step)
    for i in k:d
        step.x1[i] = step.sde.f(i,step.x0,_par(step),_t0(step))*Δt
    end
end

function eval_driftstep_xI_sym(sde::AbstractSDE{d,k,m}, method::DiscreteTimeStepping{<:RK4},x,par,Δt) where {d,k,m}
    k1 = [dyn.f[i](x, par) for i in 1:d]
    k2 = [dyn.f[i](x + Δt*k1*1//2, par) for i in 1:d]
    k3 = [dyn.f[i](x + Δt*k2*1//2, par) for i in 1:d]
    k4 = [dyn.f[i](x + Δt*k3, par) for i in 1:d]
    x + 1//6*Δt*(k1+2k2+2k3+k4)
end

function missing_states_driftstep!(step::SDEStep{d,k,m,sdeT,DiscreteTimeStepping{TDrift,TDiff}}; max_iter = 100, atol = sqrt(eps()), kwargs...) where {d,k,m,sdeT,TDrift, TDiff}
    i = 1
    while x_change < atol || i< max_iter
        iterate_xI!(step)
        x_change = norm(step.x0[j] - x for (j,x) in enumerate(step.steptracer.temp))
        update_drift_x!(step)
        i = i + 1
    end
    step
end
