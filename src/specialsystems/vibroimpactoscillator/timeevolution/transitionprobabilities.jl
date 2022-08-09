function transitionprobability(step::SDEStep{d,d,m,sdeT,method},x) where {d,m,sdeT<:SDE_VIO, method<:DiscreteTimeStepping{TDrift, TDiff}} where {TDrift, TDiff<:Maruyama}
    g = get_g(step.sde)
    r = get_wall(step)
    σ2 = _Δt0i(step) * (g(d, step.x0,_par(step),_t0(step))^ 2)
    σ2 = σ2 * (r(step.xi[2])^2)
    σ2 = σ2 + _Δti1(step) * (g(d, step.xi2,_par(step),_ti(step))^ 2) 
    # println(σ2)
    # detJ_correction = _detJ(step, x)
    normal1D_σ2(step.x1[d], σ2, x[d]);
end
