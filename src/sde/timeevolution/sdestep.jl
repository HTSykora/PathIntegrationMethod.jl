function SDEStep(sde, method, x0, x1, dt::Number; kwargs...)
    SDEStep(sde, method, x0, x1, zero(dt), dt; kwargs...)
end
function SDEStep(sde, method, x0, x1, t::AbstractVector; kwargs...)
    SDEStep(sde, method, x0, x1, similar(t,1), similar(t,1); kwargs...)
end

function SDEStep(sde::sdeT, method::methodT, x0,x1, t0, t1; precomputelevel::pclT = PreComputeNewtonStep(), kwargs...) where {sdeT<:AbstractSDE{d,k,m}, methodT <: DiscreteTimeSteppingMethod, pclT <: PreComputeLevel} where {d,k,m}
    
    _method = DiscreteTimeStepping(sde, method)
    steptracer = precomputelevel(sde,_method,x0,x1, t0, t1)
    
    SDEStep{d,k,m,sdeT,typeof(_method),typeof(steptracer),typeof(x0),typeof(x1),typeof(t0)}(sde, _method, x0, x1, t0, t1, steptracer)
end

function (pcl::Union{PreComputeNewtonStep})(sde::AbstractSDE{1,1,1}, method::DiscreteTimeStepping{<:ExplicitDriftMethod}, x0, x1, _t0, _t1)
    nothing
end
function (pcl::PreComputeNewtonStep)(sde::AbstractSDE{d,k,m}, method::DiscreteTimeStepping{<:ExplicitDriftMethod}, x0, x1, _t0, _t1) where {d,k,m}

    # @variables x[1:d] y[1:k-1] par[1:length(_par(sde))] t0 t1
    @variables x[1:d] y[1:d] par[1:length(_par(sde))] t0 t1
    
    step_sym = eval_driftstep_xI_sym(sde,method,x,par,t0,t1)
    y_eq = [y[i] - step_sym[i] for i in 1:d]
    J_sym = Symbolics.jacobian(y_eq, collect(x)); 
    JI_sym = collect(J_sym[1:k-1,1:k-1]) |> lu
    _corr = JI_sym\collect(y_eq[1:k-1]);
    x_new = [x[i] - _corr[i] for i in 1:k-1]
    detJI_inv = (det(JI_sym)^(-1))
    
    detJI⁻¹ = build_function(detJI_inv,x, par, t0, t1, expression = Val{false})
    _, xI_0! = build_function(x_new, x, collect(y[1:k-1]), par, t0, t1, expression = Val{false})
    
    _corr = lu(J_sym)\y_eq;
    x_new = [x[i] - _corr[i] for i in 1:d]
    _, x_0! = build_function(x_new, x, y, par, t0, t1, expression = Val{false})
    
    # _, xII_1! = build_function(step_sym[k:d], x, par, dt, expression = Val{false})
    SymbolicNewtonStep(xI_0!, x_0!, detJI⁻¹, similar(x0,k-1), similar(x0))
end

function (pcl::PreComputeLU)(sde::AbstractSDE{d,k,m}, method::DiscreteTimeStepping, x0, x1, Δt) where {d,k,m}
    error("LU precomputation not implemented!")
    # @variables x[1:d] y[1:k-1] par[1:length(getpar(sde))] t0 t1
    # step_sym = eval_driftstep_xI_sym(sde,method,x,par,t0,t1)
    # yI_eq = y - step_sym
    # J_sym = simplify(Symbolics.jacobian(yI_eq, x[1:k-1]), expand = true)
    
    # ...

    # StepJacobianLU(...)
end

function (pcl::PreComputeJacobian)(sde::AbstractSDE{d,k,m}, method::DiscreteTimeStepping{<:ExplicitDriftMethod}, x0, x1, Δt) where {d,k,m}
    @variables x[1:d] y[1:k-1] par[1:length(getpar(sde))] t0 t1
    step_sym = eval_driftstep_xI_sym(sde,method,x,par,t0,t1)
    yI_eq = y - step_sym
    J_sym = simplify(Symbolics.jacobian(yI_eq, x[1:k-1]), expand = true)
    _, J! = build_function(J_sym, x, par, dt, expression = Val{false})

    StepJacobian(J!, similar(x1,k-1,k-1), similar(x1,k-1))
end
# Update xI

function iterate_xI0!(step::SDEStep{d,k,m, sdeT, methodT,tracerT}) where {d, k,m, sdeT, methodT,tracerT<:SymbolicNewtonStep}
    step.steptracer.xI_0!(step.steptracer.tempI, step.x0,step.x1,_par(step),_t0(step),_t1(step))
end
function iterate_x0!(step::SDEStep{d,k,m, sdeT, methodT,tracerT}) where {d, k,m, sdeT, methodT,tracerT<:SymbolicNewtonStep}
    step.steptracer.x_0!(step.steptracer.temp, step.x0,step.x1,_par(step),_t0(step),_t1(step))
end


# Utility
_Δt(step) = _t1(step) - _t0(step)
_t0(step::SDEStep{d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT}) where {d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT} = step.t0[1] # ΔtT<:Number ???
# _t0(step::SDEStep{d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT}) where {d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT<:AbstractArray} = step.t0[1]
_t1(step::SDEStep{d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT}) where {d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT} = step.t1[1] # ΔtT<:Number ???
# _t1(step::SDEStep{d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT}) where {d, k, m, sdeT, methodT,tracerT,x0T,x1T,ΔtT<:AbstractArray} = step.t1[1]
_par(step) = step.sde.par