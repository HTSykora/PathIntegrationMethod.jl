struct TransitionTensor{N,k,NN,T,probT,pdT,tpdMX_type,dtT} <: AbstractArray{T,NN}
    sde::probT
    pdgrid::pdT
    tpdMX::tpdMX_type
    Δt::dtT
end

function TransitionTensor(sde::AbstractSDE{N,k},xs::AbstractVector,Δt::Real; kwargs...)  where {N,k}
    TransitionTensor(sde, PDGrid(sde, xs;  kwargs...), Δt; kwargs...)
end

function TransitionTensor(sde::AbstractSDE{N,k},pdgrid::PDGrid{N,k,T,xeT,xT,pT},Δt::Real) where {N,k,T,xeT,xT,pT}
    TransitionTensor{N,k,2N,eltype(pT),typeof(sde),typeof(pdgrid),Nothing,typeof(Δt)}(sde,pdgrid,nothing,Δt)
end

function TransitionTensor(sde::SDE_Oscillator1D,xs::Vector{xT},Δt::Real; kwargs...)  where {N,k,xT<:AbstractVector}
    TransitionTensor(sde, PDGrid(sde, xs;  kwargs...), Δt; kwargs...)
end
function TransitionTensor(sde::SDE_Oscillator1D,pdgrid::PDGrid{2,2,T,xeT,xT,pT},Δt::Real) where {T,xeT,xT,pT}
    TransitionTensor{2,2,2,eltype(pT),typeof(sde),typeof(pdgrid),Nothing,typeof(Δt)}(sde,pdgrid,nothing,Δt)
end


function Base.getindex(TM::TransitionTensor{N,k,NN,T,probT,pdT,tpdMX_type},idx...) where tpdMX_type<: AbstractArray{T,NN} where {N,k,NN,T,probT,pdT}
    getindex(TM.tpdMX,idx...)
end
function Base.getindex(TM::TransitionTensor{N,k,NN,T,probT,pdT,Nothing},idx...) where {N,k,NN,T,probT,pdT}
    zero(T)
end

function Base.size(TM::TransitionTensor{N,k,NN}) where {N,k,NN}
    (Base.@_inline_meta; ntuple(M -> size(TM.pdgrid.xs[mod1(M,(NN+1)÷2)], 1), Val(NN))::Dims)
end
function Base.size(TM::TransitionTensor{N,k,2}) where {N,k}
    lp = length(TM.pdgrid.p)
    (lp,lp)
end

function Base.collect(TT::TransitionTensor{N,k,NN,T,probT,pdT,Nothing,dtT}) where {N,k,NN,T,probT,pdT,dtT}
    tpdMX = zeros(T,size(TT)...);
    fillup!(tpdMX,TT)
    return TransitionTensor{N,k,NN,T,probT,pdT,typeof(tpdMX),dtT}(TT.sde, TT.pdgrid, tpdMX, TT.Δt)
end

function advance!(tt::TransitionTensor{1,1,2,T,probT,pdT,tpdMX_tpye}) where {N,k,T,probT,pdT,tpdMX_tpye<:Matrix{T}}
    mul!(tt.pdgrid.p_temp, tt.tpdMX, tt.pdgrid.p)
    tt.pdgrid.p .= sum(tt.pdgrid.p_temp)
    tt
end