#######################################
# Constructurs
function LinearInterpolation(Δ; Q_equidistant = false) 
    _Δ = Q_equidistant ? Δ : nothing
    LinearInterpolation{Q_equidistant,Nothing,typeof(_Δ)}(nothing,_Δ)
end


# Interpolation functions
function (itp::LinearInterpolation)(p::Vp,xs::Vx,x::xT, i::Integer) where {Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
    if i>=length(p)
        return p[end]
    elseif i<=1
        return p[1]
    else
        return linearinterpolate(p[i],p[i+1], (x-xs[i])/(xs[i+1] - xs[i]))
    end
end
function (itp::LinearInterpolation{true})(p::Vp,xs::Vx,x::xT, i::Integer) where {Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
    if i>=length(p)
        return p[end]
    elseif i<=1
        return p[1]
    else
        return linearinterpolate(p[i],p[i+1], (x-xs[i])/(itp.Δ))
    end
end

function linearinterpolate(p₀,p₁,θ)
    p₀ + (p₁-p₀)*θ
end

# function getweight(TT::TransitionTensor{N,k,2,T,probT,pdT},i,j) where {T,probT,pdT<:PDGrid{N,k,T2,NTuple{M,xT}}} where {T2,M,xT<:Axis{itpT}} where {itpT<:LinearInterpolation{true}} where {N,k}
#     return (j == 1 || j == length(TT.pdgrid.xs[i])) ? 0.5 : 1.0
# end
# function getarea(TT::TransitionTensor{N,k,2,T,probT,pdT}) where {T,probT,pdT<:PDGrid{N,k,T2,xsT}} where {T2,xsT<:NTuple{M,aT}} where {M,aT<: Axis{itpT}} where {itpT<:LinearInterpolation{true}} where {N,k}
#     return prod(x.itp.Δ for x in TT.pdgrid.xs) 
# end