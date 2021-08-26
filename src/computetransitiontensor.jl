function computeintegrationmatrix(sde::AbstractSDE{N,k},pdgrid,ts,method; kwargs...) where {N,k}
    warn("Pre-computation of the transitional matrices of an $(N)-dimensional problem with $(N-k+1) noises is not available yet, thus `nothing` is returned!")
    return nothing
end

function computeintegrationmatrix(sde::AbstractSDE,pdgrid::PDGrid{N,k,T},ts,method; debug_mode = false, kwargs...) where {N,k,T}
    # Prepare IK
    tpdMX = initialize_transitionmatrix(T,ts,length(pdgrid))
    
    temp = similar(pdgrid)
    IK = initialize_IK(sde,pdgrid,ts,method,temp)
    if debug_mode
        return IK
    end
    # Fill the matrix representation of the transition tensor (tpdMX)
    idx_it = get_iterator(pdgrid)
    fill_tpdMX_ts!(tpdMX,IK,idx_it,ts; kwargs...)
    tpdMX
end

# utility functions for computations
@inline initialize_transitionmatrix(T, ts::AbstractVector{eT}, l) where eT<:Number = [zeros(T,l,l) for _ in 1:(length(ts)-1)]
@inline initialize_transitionmatrix(T, ts::Number, l) = zeros(T,l,l)

@inline function initialize_IK(sde::AbstractSDE{N,N},pdgrid,ts,method,temp) where N
    t₀, t₁ = _t01(ts)
    wallID = _wallID(sde,pdgrid)
    IntegrationKernel(sde,nothing,pdgrid.xs[N],nothing,zeros(Int,N),pdgrid, t₀, t₁, method,temp, wallID); # initialize IK
end

@inline _t01(ts::Number) = zero(ts), ts
@inline _t01(ts::AbstractVector{eT}) where eT<:Number = [ts[1]], [ts[2]]
@inline _wallID(sde,pdgrid) = nothing
@inline _wallID(sde::SDE_VI_Oscillator1D,pdgrid) = ImpactInterval(zeros(eltype(pdgrid.xs[2]),2),[0],falses(1))

@inline get_iterator(pdgrid::PDGrid) = Base.Iterators.product(eachindex.(pdgrid.xs)...)
@inline get_iterator(pdgrid::PDGrid{1}) = eachindex(pdgrid.xs[1]) 

function fill_tpdMX_ts!(tpdMX,IK,idx_it,ts::Number; kwargs...)
    fill_tpdMX!(tpdMX,IK,idx_it; kwargs...)
    tpdMX
end

function fill_tpdMX_ts!(tpdMX,IK,idx_it,ts::AbstractVector{eT}; kwargs...) where {eT<:Number}
    for (jₜ) in 1:length(ts)-1
        IK.t₀[1] = ts[jₜ]
        IK.t₁[1] = ts[jₜ+1]
        fill_tpdMX!(tpdMX[jₜ],IK,idx_it; kwargs...)
    end
    tpdMX
end
function fill_tpdMX!(tpdMX,IK,idx_it; kwargs...)
    for (i,idx₁) in enumerate(idx_it)
        update_idx1!(IK,idx₁)
        get_IK_weights!(IK; integ_limits = get_integ_limits(IK), kwargs...)
        fill_to_tpdMX!(tpdMX,IK,i)
    end
    tpdMX
end
@inline function get_integ_limits(IK)
    (IK.xs[1], IK.xs[end])
end

@inline function update_idx1!(IK::IntegrationKernel{sdeT,iT0,iT1},idx₁::NTuple{M,eT}) where {sdeT,iT0,iT1<:AbstractVector{eT}} where eT<:Number where M
    for m in 1:M
        IK.idx₁[m] = idx₁[m]
    end
    nothing
end

@inline function update_idx1!(IK::IntegrationKernel{sdeT,iT0,iT1},idx₁::eT) where {sdeT<:AbstractSDE{1,1},iT0,iT1<:AbstractVector{eT}} where eT<:Number
    IK.idx₁[1] = idx₁
    nothing
end

@inline function fill_to_tpdMX!(tpdMX,IK,i)
        for j in eachindex(IK.temp)
            tpdMX[i,j] = IK.temp[j] 
        end
        nothing
end