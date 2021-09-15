function SDEStep(sde::sdeT, method::methodT, x0,x1,Δt; precomputelevel::pclT = PreComputeNewtonStep()) where {sdeT<:AbstractSDE{d,k,m}, methodT <: DiscreteTimeSteppingMethod, pclT <: PreComputeLevel} where {d,k,m}
    
    _method = DiscreteTimeStepping(method)
    steptracer = precomputelevel(sde,_method,x0,x1,Δt)
    
    SDEStep{d,k,m,sdeT,typeof(_method),typeof(steptracer),typeof(x0),typeof(x1),typeof(Δt)}(sde, _method, x0, x1, Δt, steptracer)
end

function (pcl::PreComputeNewtonStep)(sde::AbstractSDE{d,k,m}, method::DiscreteTimeSteppingMethod, x0, x1, Δt) where {d,k,m}
    @variables x[1:d] y[1:k-1] par[1:length(getpar(sde)) t0 t1
    step_sym = evalstep_sym(sde,method,x,par,t0,t1)
    yI_eq = y - step_sym[1:k-1]
    J_sym = simplify(Symbolics.jacobian(yI_eq, x[1:k-1]), expand = true)
    x_new = simplify(x[1:d] - J_sym\yI_eq, expand = true)
    detJ_inv = simplify(det(J_sym), expand = true)^(-1)
    detJ⁻¹ = build_function(detJ_inv, x, par, t0, t1, expression = Val{false})
    _, xI_0! = build_function(x_new, x, y, par, t0, t1, expression = Val{false})
    # _, xII_1! = build_function(step_sym[k:d], x, par, dt, expression = Val{false})

    SymbolicNewtonStep(xI_0!, detJ⁻¹, similar(x0,k-1))
end

function (pcl::PreComputeLU)(sde::AbstractSDE{d,k,m}, method::DiscreteTimeSteppingMethod, x0, x1, Δt) 
    error("LU precomputation not implemented!")
    # @variables x[1:d] y[1:k-1] par[1:length(getpar(sde)) t0 t1
    # step_sym = eval_driftstep_xI_sym(sde,method,x,par,t0,t1)
    # yI_eq = y - step_sym
    # J_sym = simplify(Symbolics.jacobian(yI_eq, x[1:k-1]), expand = true)
    
    # ...

    # StepJacobianLU(...)
end

function (pcl::PreComputeJacobian)(sde::AbstractSDE{d,k,m}, method::DiscreteTimeSteppingMethod, x0, x1, Δt) where {d,k,m}
    @variables x[1:d] y[1:k-1] par[1:length(getpar(sde)) t0 t1
    step_sym = eval_driftstep_xI_sym(sde,method,x,par,t0,t1)
    yI_eq = y - step_sym
    J_sym = simplify(Symbolics.jacobian(yI_eq, x[1:k-1]), expand = true)
    _, J! = build_function(J_sym, x, par, dt, expression = Val{false})

    StepJacobian(J!, similar(x1,k-1,k-1), similar(x1,k-1))
end
# Update xI

iterate_xI!(step::SDEStep{d, k, m, sdeT, methodT,tracerT}) where {d, k, m, sdeT, methodT,tracerT<:SymbolicNewtonStep}
    xI_0!(step.steptracer.temp, step.x0,step.x1,_par(step),_t0(step),_t1(step))
end

# Utility
_Δt(step) = _t1(step) - _t0(step)
_t0(step::SDEStep{d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT}) where {d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT} = step.t0[1] # ΔtT<:Number ???
# _t0(step::SDEStep{d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT}) where {d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT<:AbstractArray} = step.t0[1]
_par(step) = step.sde.par