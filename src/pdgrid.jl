# TODO: Integrate interpolation into the PDGRid!!!
Base.getindex(pdg::PDGrid, idx...) = pdg.p[idx...]
Base.size(pdg::PDGrid) = size(pdg.p)

function PDGrid(sde::AbstractSDE{N,k},_xs...;
    Q_initialize = true, axis_temp = true, kwargs...) where {N,k}
    if length(_xs) >1 && axis_temp
        xs = create_temp_axis.(Ref(Float64),_xs)
    else
        xs = _xs
    end
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
PDGrid(xs...; kwargs...) = PDGrid(SDE(),xs...; kwargs...)
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

function (p::PDGrid{1,1})(x::Number)
    p.xs[1].itp(p.p,p.xs[1],x)
end
# function (p::PDGrid{N,N})(x,y...) where N
#     idxs = [getidxs(p.xs[j+1],_y) for (j,_y) in enumerate(y)]
#     p.xs[1].itp(view(p.p,:,idxs...),p.xs[1],x)
# end

# function (p::PDGrid{N,N})(x,y...) where N
#     p.xs[1].itp(view(p.p,:,idxs...),p.xs[1],x)
# end
function (p::PDGrid)(x...)
    # TODO...
    interpolate_MV(p.p,p.xs,x...)
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
