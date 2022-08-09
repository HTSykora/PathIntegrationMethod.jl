defaultdiscreteintegrator(::SDE_VIO, di_N = 31) = Tuple(GaussLegendreIntegrator(di_N) for _ in 1:2)

Base.getindex(di::NonSmoothDiscreteIntegrator, idx...) = di.discreteintegrators[idx...]
Base.size(::NonSmoothDiscreteIntegrator{n,NoDyn}) where {n,NoDyn}= NoDyn
function DiscreteIntegrator(discreteintegrators, sdestep::NonSmoothSDEStep, res_prototype, axes::GA; kwargs...) where GA
    @assert length(discreteintegrators) == length(sdestep.sdesteps) "ERROR: length(discreteintegrators) != length(sdestep.sdesteps)"
    discreteintegrators = Tuple(DiscreteIntegrator(discreteintegrator,res_prototype, axes; kwargs...) for discreteintegrator in discreteintegrators)
    
    NonSmoothDiscreteIntegrator{1,number_of_sdesteps(sdestep),typeof(discreteintegrators)}(discreteintegrators)
end