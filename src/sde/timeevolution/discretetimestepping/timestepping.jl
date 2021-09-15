function DiscreteTimeStepping(drift)
    DiscreteTimeStepping(drift,Maruyama())
end

DiscreteTimeStepping(method::DiscreteTimeStepping) = method



# Evaluate time stepping


# Utils
function (J::StepJacobian)(args::Vararg{Any,N}) where {N}
    J.J!(J.JM,args...)
    J.JM
end