function PathIntegrationProblem(sde::AbstractSDE{N,k},pdgrid::PDGrid{N,k,T,xT,pT},Δt::Real; method = EulerMaruyama(), precompute = false, kwargs...) where {N,k,T,xT,pT}
    if precompute
        error("Pre-computation of the transitional matrices of an $(N)-dimensional problem with $(N-k+1) noises is not available yet")
    else
        return PathIntegrationProblem{N,k,typeof(sde),typeof(pdgrid),Nothing,typeof(Δt),typeof(method)}(sde,pdgrid,nothing,Δt, method)
    end
end

function PathIntegrationProblem(sde::AbstractSDE{N,k},Δt::Real,xs...; kwargs...)  where {N,k}
    PathIntegrationProblem(sde, PDGrid(sde, xs...;  kwargs...), Δt; kwargs...)
end


# function PathIntegrationProblem(sde::SDE_Oscillator1D,pdgrid::PDGrid{2,2},Δt::Real; method = EulerMaruyama())
#     PathIntegrationProblem{2,2,typeof(sde),typeof(pdgrid),Nothing,typeof(Δt),typeof(method)}(sde,pdgrid,nothing,Δt, method)
# end