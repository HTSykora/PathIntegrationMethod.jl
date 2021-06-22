struct TransitionTensor{T,N,probT,pdT,dtT} <: AbstractArray{T,N}
    sde::probT
    pdgrid::pdT
    Δt::dtT
end

# TPDMX::TPDMX_type
# Δt::Δt_type

# function _tp(prob::SDEProblem,x,t,x0,t0) # transition probability
#     σ² = prob.g(x0)^2*(t-t0)
#     μ = x0 + prob.f(x0)*(t-t0)
#     exp(-0.5*((x - μ)^2 / σ²))/sqrt(2π*σ²)
# end


# function TransitionTensor(prob::SDEProblem,pdg::pdgT,Δt) where pdgT<:PDGrid{T,xT} where xT<: AbstractVector{T} where T<:Number
#     TransitionTensor{T,2,typeof(prob),pdgT,typeof(Δt)}(prob,pdg,Δt)
# end
# function TransitionTensor(prob::SDEProblem,x::AbstractVector{T},Δt) where T <: Number
#     p = [zero(x)]; 
#     lx = length(x);
#     if isodd(lx)
#         for i in 1:3
#             p[1][lx÷2 + i] = T(1/3);
#         end
#     else
#         for i in 1:2
#             p[1][lx÷2 + i] = T(0.5)
#         end
#     end
#     pdg = PDGrid(x[2]-x[1],x,p);
#     TransitionTensor{T,2,typeof(prob),typeof(pdg),typeof(Δt)}(prob,pdg,Δt)
# end

# function Base.size(t::TransitionTensor{T,2}) where T
#     l = length(t.pdg.x)
#     (l,l)
# end
# function Base.getindex(t::TransitionTensor{T,2},i,j) where T
#     _tp(t.prob,t.pdg.x[i],t.Δt,t.pdg.x[j],0.)*t.pdg.Δx # transition probability
# end