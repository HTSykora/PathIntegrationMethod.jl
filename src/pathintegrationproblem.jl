function PathIntegrationProblem(sde::AbstractSDE{N,k}, pdgrid::PDGrid{N,k,T,xT,pT}, ts; method = EulerMaruyama(), precompute = false, kwargs...) where {N,k,T,xT,pT}
    # ts == Δt if ts::Number
    # ts: i-th time interval between [ ts[i], ts[i+1] ] if ts::{Tuple/AbstractVector}
    if ts isa Number
        step_idx = nothing
    else
        step_idx = [0]
    end

    if precompute
        tpdMX = computeintegrationmatrix(sde,pdgrid, ts, method; kwargs...)
    else
        tpdMX = nothing 
    end
    pp_ts = _ts(ts)
    return PathIntegrationProblem{N,k,typeof(sde),typeof(pdgrid),typeof(tpdMX),typeof(pp_ts),typeof(method),typeof(step_idx)}(sde, pdgrid, tpdMX, pp_ts, method, step_idx)
end

function PathIntegrationProblem(sde::AbstractSDE{N,k},ts, xs...; kwargs...)  where {N,k}
    PathIntegrationProblem(sde, PDGrid(sde, xs...;  kwargs...), ts; kwargs...)
end

_ts(ts) = ts
_ts(ts::Number) = (zero(ts), ts)

##############################################
# Advance functions
function advance!(pip::PathIntegrationProblem{N,k,sdeT,pdT,tpdMX_tpye}) where {N,k,sdeT,pdT,tpdMX_tpye<:mT} where mT<: AbstractMatrix{T} where T<:Number
    mul!(vec(pip.pdgrid.p_temp), pip.tpdMX, vec(pip.pdgrid.p))
    pip.pdgrid.p .= pip.pdgrid.p_temp ./ _integrate(pip.pdgrid.p_temp, pip.pdgrid.xs)
    pip
end
function advance!(pip::PathIntegrationProblem{N,k,sdeT,pdT,tpdMX_tpye}) where {N,k,sdeT,pdT,tpdMX_tpye<:vmT} where vmT<: AbstractVector{mT} where mT<:AbstractMatrix{T} where T<:Number
    pip.step_idx[1] += mod1(pip.step_idx[1] + 1, length(pip.tpdMX))
    mul!(vec(pip.pdgrid.p_temp), pip.tpdMX[pip.step_idx[1]], vec(pip.pdgrid.p))
    pip.pdgrid.p .= pip.pdgrid.p_temp ./ _integrate(pip.pdgrid.p_temp, pip.pdgrid.xs)
    pip
end

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
