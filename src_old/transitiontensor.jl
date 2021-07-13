function Base.getindex(TM::PathIntegrationProblem{N,k,NN,T,probT,pdT,tpdMX_type},idx...) where tpdMX_type<: AbstractArray{T,NN} where {N,k,NN,T,probT,pdT}
    getindex(TM.tpdMX,idx...)
end
function Base.getindex(TM::PathIntegrationProblem{N,k,NN,T,probT,pdT,Nothing},idx...) where {N,k,NN,T,probT,pdT}
    zero(T)
end






function PathIntegrationProblem(sde::SDE_Oscillator1D,Δt::Real,xs...; kwargs...)  where {N,k}
    PathIntegrationProblem(sde, PDGrid(sde, xs...;  kwargs...), Δt; kwargs...)
end




function Base.size(TM::PathIntegrationProblem{N,k,NN}) where {N,k,NN}
    (Base.@_inline_meta; ntuple(M -> size(TM.pdgrid.xs[mod1(M,(NN+1)÷2)], 1), Val(NN))::Dims)
end
function Base.size(TM::PathIntegrationProblem{N,k,2}) where {N,k}
    lp = length(TM.pdgrid.p)
    (lp,lp)
end

function Base.collect(TT::PathIntegrationProblem{N,k,NN,T,probT,pdT,Nothing,dtT,methodT}) where {N,k,NN,T,probT,pdT,dtT,methodT}
    tpdMX = zeros(T,size(TT)...);
    fillup!(tpdMX,TT)
    return PathIntegrationProblem{N,k,NN,T,probT,pdT,typeof(tpdMX),dtT,methodT}(TT.sde, TT.pdgrid, tpdMX, TT.Δt,TT.method)
end

function advance!(tt::PathIntegrationProblem{1,1,2,T,probT,pdT,Nothing}) where {T,probT,pdT}
    zΔt = zero(tt.Δt)
    for (i,x₁) in enumerate(tt.pdgrid.xs[1])
        tt.pdgrid.p_temp[i] = _integrate(x₀-> _tp(tt.sde,x₁,tt.Δt, x₀,zΔt, method = tt.method),tt.pdgrid.p,tt.pdgrid.xs)
    end
    tt.pdgrid.p .= tt.pdgrid.p_temp ./ _integrate(tt.pdgrid.p_temp, tt.pdgrid.xs)
    # tt.pdgrid.p .= tt.pdgrid.p_temp
    tt
end
function get_tp(tt::PathIntegrationProblem{1,1,2,T,probT,pdT,Nothing},x₁) where {T,probT,pdT}
    zΔt = zero(tt.Δt)
    [_tp(tt.sde,x₁,tt.Δt, x₀,zΔt, method = tt.method) for x₀ in tt.pdgrid.xs[1]]
    # [_tp(tt.sde,x₀,tt.Δt, x₁,zΔt, method = tt.method) for x₀ in tt.pdgrid.xs[1]]
end

function advance!(tt::PathIntegrationProblem{1,1,2,T,probT,pdT,tpdMX_tpye}) where {N,k,T,probT,pdT,tpdMX_tpye<:Matrix{T}}
    mul!(tt.pdgrid.p_temp, tt.tpdMX, tt.pdgrid.p)
    tt.pdgrid.p .= tt.pdgrid.p_temp ./ _integrate(tt.pdgrid.p_temp, tt.pdgrid.xs[1])
    tt
end

function advance!(tt::PathIntegrationProblem{2,2,2,T,probT,pdT,tpdMX_tpye}) where {N,k,T,probT,pdT,tpdMX_tpye<:Matrix{T}}
    for (j,v) in enumerate(tt.pdgrid.xs[2])
        for (i,x) in enumerate(tt.pdgrid.xs[1])
            tt.pdgrid.p_temp[i,j] = tt.pdgrid(tt.pdgrid.ξ_temp[i,j],v)
            # tt.pdgrid.p_temp[i,j] = tt.pdgrid.p[i,j]
        end
    end
    # mul!(vec(tt.pdgrid.p), tt.tpdMX, vec(tt.pdgrid.p_temp))
    vec(tt.pdgrid.p).= tt.tpdMX*vec(tt.pdgrid.p_temp)
    tt.pdgrid.p ./= sum(tt.pdgrid.p) # TODO
    tt
end

# function advance!(tt::TransitionTensor{2,2,2,T,probT,pdT,Nothing}) where {T,probT,pdT}
#     for (j,v) in enumerate(tt.pdgrid.xs[2])
#         for (i,x) in enumerate(tt.pdgrid.xs[1])
#             tt.pdgrid.p_temp[i,j] = quadgk(v₀ -> _TP(tt.sde,(x,v),tt.Δt,(get_ξ(tt.method,tt.sde,tt.Δt,0.,x,v₀),v₀),0.)*tt.pdgrid(get_ξ(tt.method,tt.sde,tt.Δt,0.,x,v₀),v₀),tt.pdgrid.xs[2][1],tt.pdgrid.xs[2][end])[1]
#         end
#     end
#     tt.pdgrid.p .= tt.pdgrid.p_temp ./ _integrate(tt.pdgrid.p_temp,tt.pdgrid.xs)
# end
function advance!(tt::PathIntegrationProblem{2,2,2,T,probT,pdT,Nothing}) where {T,probT,pdT}
    for (j,v) in enumerate(tt.pdgrid.xs[2])
        for (i,x) in enumerate(tt.pdgrid.xs[1])
            for (j₀,v₀) in enumerate(tt.pdgrid.xs[2])
                tt.pdgrid.i_temp[j₀] = _TP(tt.sde,(x,v),tt.Δt,(get_ξ(tt.method,tt.sde,tt.Δt,0.,x,v₀),v₀),0.)*tt.pdgrid(get_ξ(tt.method,tt.sde,tt.Δt,0.,x,v₀),v₀)
            end
            tt.pdgrid.p_temp[i,j] = _integrate(tt.pdgrid.i_temp,tt.pdgrid.xs[2])
        end
    end
    tt.pdgrid.p .= tt.pdgrid.p_temp ./ _integrate(tt.pdgrid.p_temp,tt.pdgrid.xs)
end