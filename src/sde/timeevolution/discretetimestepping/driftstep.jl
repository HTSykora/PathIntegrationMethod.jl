function eval_driftstep!(step::SDEStep{d,k,m, sdeT, methodT}) where {d,k,m,sdeT, methodT<: DiscreteTimeStepping{TDrift}} where {TDrift<:Euler}
    Δt = _Δt(step)
    for i in 1:d
        step.x1[i] = step.x0[i] + step.sde.f(i,step.x0,_par(step),_t0(step))*Δt
    end
end
function update_x1_kd!(step::SDEStep{d,k,m, sdeT, methodT}) where {d,k,m,sdeT, methodT<: DiscreteTimeStepping{TDrift}} where {TDrift<:Euler}
    Δt = _Δt(step)
    for i in k:d
        step.x1[i] = step.x0[i] + step.sde.f(i,step.x0,_par(step),_t0(step))*Δt
    end
end

function eval_driftstep!(step::SDEStep{d,k,m, sdeT, methodT}) where {d,k,m,sdeT, methodT<: DiscreteTimeStepping{TDrift}} where {TDrift<:RungeKutta{ord}} where ord
    _eval_driftstep!(step)
    fill_to_x1(step)
end
function _eval_driftstep!(step,x0,x1)
    Δt = _Δt(step)
    for i in 1:d
        step.method.ks[1][i] = x0[i] + step.sde.f(i,x0,_par(step),_t0(step))*Δt
    end

    for j in 1:ord-1
        step.method.temp .= step.method.x0
        for a in step.method.BT.a[j]
            step.method.temp .= step.method.temp .+ a.weight .* a.val
        end
        tj = _t0(step)+Δt*(1+step.method.BT.c[j])
        for i in 1:d
            step.method.ks[j+1][i] = x0[i] + step.sde.f(i,step.method.temp,_par(step),tj)*Δt
        end
    end
end
function fill_to_x1(step)
    step.x1 .= step.x0
    for b in step.method.BT.b
        step.x1 .= step.x1 .+ b.weight .* b.val
    end
end
function fill_to_x1(step::SDEStep{d,k,m, sdeT, methodT},i) where {d,k,m,sdeT, methodT<: DiscreteTimeStepping{TDrift}} where {TDrift<:RungeKutta{ord}} where ord
    step.x1 .= step.x0
    for b in step.method.BT.b
        for j in i:d
            step.x1[j] .= step.x1[j] .+ b.weight .* b.val[j]
        end
    end
end
function update_x1_kd!(step::SDEStep{d,k,m, sdeT, methodT}) where {d,k,m,sdeT, methodT<: DiscreteTimeStepping{TDrift}} where {TDrift<:RungeKutta{ord}} where ord
    _eval_driftstep!(step,step.x0,step.x1_temp)
    for i in k:d
        step.x1[i] = step.x1_temp[i]
    end
end


function eval_driftstep_xI_sym(sde::AbstractSDE{d,k,m}, method::DiscreteTimeStepping{<:Euler}, x, par, t0, t1) where {d,k,m}
    [x[i] + sde.f(i,x,par,t0)*(t1-t0) for i in 1:d]
end
function eval_driftstep_xI_sym(sde::AbstractSDE{d,k,m}, method::DiscreteTimeStepping{<:RungeKutta{ord}},x,par,t0,t1) where {d,k,m,ord}
    Δt = t1 - t0
    ks = [[x[i] + sde.f(i,x,par,t0)*(t1-t0) for i in 1:d]]
    
    temp = similar(x)
    for j in 1:ord-1
        temp .= x
        for a in method.BT.a[j]
            temp .= temp .+ a._weight .* ks[a.idx]
        end
        tj = t0 + Δt*(1 + step.method.BT._c[j])
        push!(ks, [x[i] + sde.f(i,method.temp,par,tj)*Δt for i in 1:d])
    end
    
    temp .= x
    for b in step.method.BT.b
        temp .= temp .+ b._weight .* ks[b.idx]
    end
    
    return temp
end

function update_drift_x!(step::SDEStep{d,k,m, sdeT, methodT}) where {d,k,m,sdeT,methodT<: DiscreteTimeStepping{TDrift}} where {TDrift}
    for i in 1:k-1
        step.x0[i] = step.steptracer.temp[i]
    end
    
    update_x1_kd!(step)
end

function compute_missing_states_driftstep!(step::SDEStep{d,1,m,sdeT,DiscreteTimeStepping{TDrift,TDiff}}; kwargs...) where {d,m,sdeT,TDrift, TDiff}
    eval_driftstep!(step)
end

function compute_missing_states_driftstep!(step::SDEStep{d,k,m,sdeT,DiscreteTimeStepping{TDrift,TDiff}}; max_iter = 100, atol = sqrt(eps()), kwargs...) where {d,k,m,sdeT,TDrift, TDiff}
    i = 1
    x_change = 2atol
    while x_change > atol && i < max_iter
        iterate_xI!(step)
        x_change = norm(step.x0[j] - x for (j,x) in enumerate(step.steptracer.temp))
        update_drift_x!(step)
        i = i + 1
    end
    # println("iterations: $(i-1)")
    nothing
end

# Butcher Tableus for Runge Kutta method
BTElement(T::DataType,idx,_weight,val) = BTElement(idx,_weight,T(weight),val)
function RungeKutta(order::Integer, BT::btT,ks::ksT, temp::tT) where{btT,ksT,tT}
    RungeKutta{order,btT,ksT,tT}(BT,ks,temp)
end
function RK4(T::DataType = Float64)
    temp = Vector{T}(undef,0)
    ks = Tuple(similar(temp) for _ in 1:4)
    _c = (1//2, 1//2, 1)
    c = T.(_c)
    b = (BTElement(T,1, 1//3, ks[1]), BTElement(T,2, 1//6, ks[2]), BTElement(T,3, 1//6, ks[3]), BTElement(T,4, 1//3, ks[4]))
    a = ((BTElement(T,1,1//2, ks[1]),),(BTElement(T,2,1//2, ks[2]),),(BTElement(T,3,1, ks[3]),))
    
    RungeKutta(4, ButcherTableau(a,b,c,_c), ks, temp)
end
