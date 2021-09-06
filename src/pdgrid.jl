# TODO: Integrate interpolation into the PDGRid!!!
Base.getindex(pdg::PDGrid, idx...) = pdg.p[idx...]
Base.size(pdg::PDGrid) = size(pdg.p)
# function PDgrid(p,_xs...; axis_temp = true, kwargs...)
#     if length(_xs) >1 && axis_temp
#         xs = create_temp_axis.(Ref(Float64),_xs)
#     else
#         xs = _xs
#     end
#     return PDGrid{N,k,eltype(p),typeof(xs),typeof(p),typeof(ξ_temp),typeof(grid),typeof(i_temp)}(xs, p, similar(p), ξ_temp, grid,i_temp)
# end

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
        initialize!(p,xs...; kwargs...)
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
PDGrid(xs...; kwargs...) = PDGrid(DummySDE(length(xs),length(xs)),xs...; kwargs...)
function PDGrid(pM::AbstractArray{T,N}, xs...; kwargs...) where {T,N}
    pdgrid = PDGrid(DummySDE(N,N),xs...; Q_initialize = false, kwargs...)
    pdgrid.p .= pM
    pdgrid
end
# function PDGrid(sde::SDE{1,1},xs::xsT; kwargs...) where xsT<:Axis
#     PDGrid(sde, xs; kwargs...)
# end

@inline function initialize!(p::AbstractArray{T,N}; kwargs...) where {T<:Number,N}
    p_size=size(p)
    Q_oddp = [isodd(l) for l in p_size]
    weight = 1/prod(2 + oddp for oddp in Q_oddp)
    middims = ntuple(i->(p_size[i]-1)÷2:((p_size[i]+1)÷2 + Q_oddp[i]),Val(N))
    p[middims...] .= weight
    p
end

@inline function initialize!(p::AbstractArray{T,N},xs...; μ_init = nothing, σ_init = nothing,    kwargs...) where {T<:Number,N}
    if μ_init isa Nothing
        μs = [(x[end]+x[1])/2 for x in xs]
    elseif μ_init isa Number
        μs = [μ_init for _ in xs]
    elseif μ_init isa Vector
        @assert length(μ_init) == length(xs) "Wrong number of initial μ values are given"
        μs = μ_init
    end

    if σ_init isa Nothing
        σ²s = [((x[end]-x[1])/12)^2 for x in xs]
    elseif σ_init isa Number
        σ²s = [σ_init^2 for _ in xs]
    elseif σ_init isa Vector
        @assert length(σ_init) == length(xs) "Wrong number of initial σ values are given"
        σ²s = σ_init.^2
    end
    
    idx_it = Base.Iterators.product(eachindex.(xs)...)
    
    for idxs in idx_it
        p[idxs...] = prod(normal1D_σ2(μs[i],σ²s[i],xs[i][idx]) for (i,idx) in enumerate(idxs))
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
