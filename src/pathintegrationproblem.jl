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
function advance!(tt::PathIntegrationProblem{N,k,sdeT,pdT,tpdMX_tpye}) where {N,k,sdeT,pdT,tpdMX_tpye<:Matrix{T}} where T<:Number
    mul!(vec(tt.pdgrid.p_temp), tt.tpdMX, vec(tt.pdgrid.p))
    tt.pdgrid.p .= tt.pdgrid.p_temp ./ _integrate(tt.pdgrid.p_temp, tt.pdgrid.xs)
    tt
end
function advance!(tt::PathIntegrationProblem{2,2,sdeT,pdT,tpdMX_tpye}) where {sdeT<:SDE_Oscillator1D,pdT,tpdMX_tpye} where T<:Number
    pszs = size(tt.pdgrid)
    # for j in 1:pszs[2], i in 1:pszs[1]
    #         tt.pdgrid.p_temp[i,j] = sum(tt.pdgrid.p[k,l]*tt.tpdMX[k,l,i,j] for l in 1:pszs[2], k in 1:pszs[1])
    # end
    mul!(vec(tt.pdgrid.p_temp), tt.tpdMX, vec(tt.pdgrid.p))
    tt.pdgrid.p .= tt.pdgrid.p_temp ./ _integrate(tt.pdgrid.p_temp, tt.pdgrid.xs)
    tt
end

# function PathIntegrationProblem(sde::SDE_Oscillator1D,pdgrid::PDGrid{2,2},Δt::Real; method = EulerMaruyama())
#     PathIntegrationProblem{2,2,typeof(sde),typeof(pdgrid),Nothing,typeof(Δt),typeof(method)}(sde,pdgrid,nothing,Δt, method)
# end

function get_stationary_by_eigenvectors(pip::PathIntegrationProblem; ev_id = 1, kwargs...)
    pip_res = deepcopy(pip)
    esys = eigs(pip_res.tpdMX; kwargs);
    e_vec = esys[2][:,ev_id] .|> real
    e_val = esys[1][ev_id] .|> real
    for (i,_ev) in enumerate(e_vec)
        pip_res.pdgrid.p[i] .= _ev
    end
    pip_res.pdgrid.p .= pip_res.pdgrid.p ./ _integrate(pip_res.pdgrid)

    pip_res.pdgrid, e_val
end
