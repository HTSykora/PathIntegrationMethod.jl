struct PDGrid{N,k,T,xT,pT,ξT} <: AbstractArray{T,N}
    xs::xT
    p::pT
    p_temp::pT
    ξ_temp::ξT
end

Base.getindex(pdg::PDGrid, idx...) = pdg.p[idx...]
Base.size(pdg::PDGrid) = size(pdg.p)

# TODO: Integrate interpolation into the PDGRid!!!
function PDGrid(sde::AbstractSDE{N,k},xs::xsT; Q_initialize = true, kwargs...) where xsT<:AbstractVector{xT} where xT<:AbstractVector where {N,k}
    lens = length.(xs);
    p = zeros(eltype(xs[1]),lens...)
    if Q_initialize
        initialize!(p,xs)
    end
    if N == k
        ξ_temp = similar(p)
    end
    PDGrid{N,k,eltype(p),xsT,typeof(p),typeof(p)}(xs,p,similar(p), ξ_temp)
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
@inline function initialize!(p::AbstractArray{T,N},xs::Txs) where {T<:Number,N} where Txs<:AbstractVector{Tx} where Tx<:AbstractVector
    μs = [(x[end]+x[1])/2 for x in xs]
    σ²s = [((x[end]-x[1])/12)^2 for x in xs]
    idx_it = Base.Iterators.product(eachindex.(xs)...)
    
    for idxs in idx_it
        p[idxs...] = prod(normal1D(μs[i],σ²s[i],xs[i][idx]) for (i,idx) in enumerate(idxs))
    end
end

function (p::PDGrid{1,1})(x)
    p.xs[1].itp(p.p,p.xs[1],x)
end
function (p::PDGrid{N,N})(x,y...) where N
    idxs = [getidxs(p.xs[j+1],_y) for (j,_y) in enumerate(y)]
    p.xs[1].itp(view(p.p,:,idxs...),p.xs[1],x)
end
function getidxs(xs,x::xT) where xT
    if x isa Union{Rational,Integer}
        @inbounds i = searchsortedlast(xs, x)  
    else
        @inbounds i = searchsortedlast(xs, x - 10eps(xT))
    end
    if i>=length(xs)  # it can happen if t ≈ t0
        return length(xs)
    else
        return i+1
    end
end


## Normalize functions
function renormalize!(p::PDGrid{1,1})
    renormalize!(p.p,p.xs[1])
end
function renormalize!(p::Vp,xs::Axis) where Vp<:AbstractVector{Tp} where Tp<:Number
    p ./= integrate_p(p,xs)
end
