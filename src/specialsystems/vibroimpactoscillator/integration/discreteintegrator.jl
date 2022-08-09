Base.getindex(di::NonSmoothDiscreteIntegrator, idx...) = di.discreteintegrators[idx...]
Base.size(::NonSmoothDiscreteIntegrator{n,NoDyn}) where {n,NoDyn}= NoDyn
function DiscreteIntegrator(discreteintegrator, sdestep::NonSmoothSDEStep, res_prototype, N::Union{NTuple{1,<:Integer},<:Integer,AbstractArray{<:Integer}}, axes::GA; kwargs...) where GA
    discreteintegrators = Tuple(DiscreteIntegrator(discreteintegrator,res_prototype, N, axes; kwargs...) for _ in sdestep.sdesteps)
    
    NonSmoothDiscreteIntegrator{1,number_of_sdesteps(sdestep),typeof(discreteintegrators)}(discreteintegrators)
end