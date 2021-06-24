struct PDGrid{N,k,T,xeT,xT,pT} <: AbstractArray{T,N}
    Δx::xeT
    xs::xT
    p::pT
    p_temp::pT
    itp::itpT # Chebyshev or equidistant linear grid
end

Base.getindex(pdg::PDGrid, idx...) = pdg.p[idx...]
Base.size(pdg::PDGrid) = size(pdg.p)

function PDGrid(sde::SDE{N,k},xs::xsT; Q_equidistant = true, Q_initialize = true, kwargs...) where xsT<:AbstractVector{xT} where xT<:AbstractVector where {N,k}
    if Q_equidistant
        Δx = [x[2]-x[1] for x in xs]
    else
        Δx = nothing
    end
    lens = length.(xs);
    p = zeros(eltype(xs[1]),lens...)
    if Q_initialize
        initialize!(p)
    end
    PDGrid{N,k,eltype(p),typeof(Δx),xsT,typeof(p)}(Δx,xs,p,similar(p))
end

function PDGrid(sde::SDE{1,1},xs::xsT; kwargs...) where xsT<:AbstractVector{xT} where xT<:Number
    PDGrid(sde,[xs]; kwargs...)
end

@inline function initialize!(p::AbstractArray{T,N}) where {T<:Number,N}
    p_size=size(p)
    Q_oddp = [isodd(l) for l in p_size]
    weight = 1/prod(2 + oddp for oddp in Q_oddp)
    middims = ntuple(i->(p_size[i]-1)÷2:((p_size[i]+1)÷2 + Q_oddp[i]),Val(N))
    p[middims...] .= weight
    p
end

