abstract type AbstractAxis{T} <:AbstractVector{T} end 

struct Axis{itpT,T,xT} <: AbstractAxis{T}
    itp::itpT
    x::xT
end

Base.getindex(a::Axis,idx...) = a.x[idx...]
Base.size(a::Axis) = size(a.x)

function Axis(start,stop,num::Int; lvl = 1, interpolation = :linear)
    if interpolation == :linear
        xs = LinRange(start,stop,num)
        itp = LinearInterpolation(num,(stop-start)/(num-1),Q_equidistant = true);
    elseif interpolation == :chebyshev
        xs = chebygrid(start,stop,num)
        itp = ChebyshevInterpolation(num, start, stop);
    else
        error("Interpolation should be `:linear` or `:chebyshev`")
    end
    return Axis{typeof(itp),eltype(xs),typeof(xs)}(itp,xs)
end

struct PDGrid{N,k,T,xT,pT,ξT,gridT,iT} <: AbstractArray{T,N}
    xs::xT
    p::pT
    p_temp::pT
    ξ_temp::ξT
    grid::gridT
    i_temp::iT
end

Base.getindex(pdg::PDGrid, idx...) = pdg.p[idx...]
Base.size(pdg::PDGrid) = size(pdg.p)

# TODO: Integrate interpolation into the PDGRid!!!

function PDGrid(sde::AbstractSDE{N,k},xs...; Q_initialize = true, kwargs...) where {N,k}
    lens = length.(xs);
    p = zeros(eltype(xs[1]),lens...)
    if Q_initialize
        initialize!(p,xs...)
    end

    if N == k
        # ξ_temp = similar(p)
        ξ_temp = nothing
    else
        ξ_temp = nothing
    end

    if N>1
        i_temp = similar(p,lens[end-(N-k):end]...)
    else
        i_temp = nothing
    end
    if xs isa NTuple{N,Axis{<:LinearInterpolation{true}}} where N
        grid = RectangleGrid(xs...)
    else
        grid = nothing
    end
    return PDGrid{N,k,eltype(p),typeof(xs),typeof(p),typeof(ξ_temp),typeof(grid),typeof(i_temp)}(xs, p, similar(p), ξ_temp, grid,i_temp)
end

PDGrid(sde::AbstractSDE,xs::AV; kwargs...) where AV<:AbstractVector{AX} where AX<:Axis = PDGrid(sde,xs...; kwargs...)

# function PDGrid(sde::SDE{1,1},xs::xsT; kwargs...) where xsT<:Axis
#     PDGrid(sde, xs; kwargs...)
# end

@inline function initialize!(p::AbstractArray{T,N}) where {T<:Number,N}
    p_size=size(p)
    Q_oddp = [isodd(l) for l in p_size]
    weight = 1/prod(2 + oddp for oddp in Q_oddp)
    middims = ntuple(i->(p_size[i]-1)÷2:((p_size[i]+1)÷2 + Q_oddp[i]),Val(N))
    p[middims...] .= weight
    p
end

@inline function initialize!(p::AbstractArray{T,N},xs...) where {T<:Number,N}
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
# function (p::PDGrid{N,N})(x,y...) where N
#     idxs = [getidxs(p.xs[j+1],_y) for (j,_y) in enumerate(y)]
#     p.xs[1].itp(view(p.p,:,idxs...),p.xs[1],x)
# end

# function (p::PDGrid{N,N})(x,y...) where N
#     p.xs[1].itp(view(p.p,:,idxs...),p.xs[1],x)
# end
function (p::PDGrid{N, k, T, NTuple{M,aT}} where {N,k,T,M,aT<:Axis{<:LinearInterpolation{true}}})(x...)
    interpolate(p.grid,p.p,[x...])
end


function interpolate_1(p,x,y...)
    idxs = [getidxs(p.xs[j+1],_y) for (j,_y) in enumerate(y)]
    p.xs[1].itp(view(p.p,:,idxs...),p.xs[1],x)
end

# function (p::PDGrid{N,k,T,xT,pT,ξT,gridT},x...) where {N,k,T,xT,pT,ξT,gridT<:RectangleGrid}
#     interpolate(pdgrid.grid,pdgrid.p,x)
# end

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
