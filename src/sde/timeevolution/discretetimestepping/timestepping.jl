function DiscreteTimeStepping(drift::DiscreteTimeSteppingMethod)
    DiscreteTimeStepping(drift,Maruyama())
end

DiscreteTimeStepping(method::DiscreteTimeStepping) = method

function DiscreteTimeStepping(sde::AbstractSDE, drift::DiscreteTimeSteppingMethod)
    _method = DiscreteTimeStepping(drift,Maruyama())
    DiscreteTimeStepping(sde, _method)
end

function DiscreteTimeStepping(sde::AbstractSDE, method::DiscreteTimeStepping)
    resize_tempstates!(sde,method)
    method
end
resize_tempstates!(sde,method) = nothing
function resize_tempstates!(sde::AbstractSDE{d,k,m},method::DiscreteTimeStepping{RK}) where {d,k,m,RK<:RungeKutta}
    resize!.(method.drift.ks,d)
    resize!(method.drift.temp,d)
    method
end

# Evaluate time stepping


# Utils
function (J::StepJacobian)(args::Vararg{Any,N}) where {N}
    J.J!(J.JM,args...)
    J.JM
end