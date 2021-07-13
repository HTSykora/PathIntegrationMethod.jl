function PathIntegrationProblem(sde::AbstractSDE{N,k},pdgrid::PDGrid{N,k,T,xT,pT},Δt::Real; method = EulerMaruyama(), precompute = false, kwargs...) where {N,k,T,xT,pT}
    if precompute
        tpdMX = computeintegrationmatrix(sde,pdgrid,Δt,method; kwargs...)
    else
        tpdMX = nothing
    end
    return PathIntegrationProblem{N,k,typeof(sde),typeof(pdgrid),typeof(tpdMX),typeof(Δt),typeof(method)}(sde, pdgrid, tpdMX, Δt, method)
end

function PathIntegrationProblem(sde::AbstractSDE{N,k},Δt::Real,xs...; kwargs...)  where {N,k}
    PathIntegrationProblem(sde, PDGrid(sde, xs...;  kwargs...), Δt; kwargs...)
end

##############################################
# Advance functions
function advance!(tt::PathIntegrationProblem{1,1,sdeT,pdT,tpdMX_tpye}) where {sdeT,pdT,tpdMX_tpye<:Matrix{T}} where T<:Number
    mul!(tt.pdgrid.p_temp, tt.tpdMX, tt.pdgrid.p)
    tt.pdgrid.p .= tt.pdgrid.p_temp ./ _integrate(tt.pdgrid.p_temp, tt.pdgrid.xs[1])
    tt
end

# function PathIntegrationProblem(sde::SDE_Oscillator1D,pdgrid::PDGrid{2,2},Δt::Real; method = EulerMaruyama())
#     PathIntegrationProblem{2,2,typeof(sde),typeof(pdgrid),Nothing,typeof(Δt),typeof(method)}(sde,pdgrid,nothing,Δt, method)
# end