# Trapezoidal weights for linear interpolation
Base.getindex(t::TrapezoidalWeights{T},idx::Integer) where T = (idx == 1 || idx == t.l) ? T(0.5) : T(1.)
Base.size(t::TrapezoidalWeights) = (t.l,)

#######################################
# Constructurs
LinearInterpolation(Δ, n::Integer = 0; kwargs...) = LinearInterpolation(Float64,Δ,n;kwargs...)
function LinearInterpolation(T::DataType,Δ::Number, n::Integer = 0; Q_equidistant = false, create_temp = false) 
    _Δ = Q_equidistant ? Δ : nothing
    temp = create_temp ? zeros(T,n) : nothing
    LinearInterpolation{Q_equidistant,typeof(temp),typeof(_Δ)}(temp,_Δ)
end

# Reconstructor
function new_itp(T::DataType,itp::LinearInterpolation{TF,Ttemp,TΔ}, n ) where {TF,Ttemp,TΔ} # Ttemp <: Nothing
    LinearInterpolation(T,itp.Δ, n, Q_equidistant = TF, create_temp = true)
end

#######################################
# Interpolation functions
getΔ(itp::LinearInterpolation{true},xs,i) = itp.Δ
getΔ(itp::LinearInterpolation{false},xs,i) = xs[i+1] - xs[i]

function (itp::LinearInterpolation)(p::Vp,xs::Vx,x::xT, i::Integer) where {Vp<:AbstractVector{Tp},Vx<:AbstractVector{Tx},xT<:Number} where Tp<:Number where Tx<:Number
    if i>=length(p)
        return p[end]
    elseif i<=1
        return p[1]
    else
        Δ = getΔ(itp,xs,i)
        return linearinterpolate(p[i],p[i+1], (x-xs[i])/Δ)
    end
end

function linearinterpolate(p₀,p₁,θ)
    p₀ + (p₁-p₀)*θ
end

function basefun_vals!(itp::LinearInterpolation,vals,xs::Vx,x) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    # i = 0... N
    vals .= zero(eltype(vals));
    needsinterpolation, i = find_idx(xs,x)
    if needsinterpolation
        Δ = getΔ(itp,xs,i)
        θ = eltype(vals)(x-xs[i]/Δ)
        vals[i] = 1-θ
        vals[i+1] = θ
    else
        vals[i] = one(eltype(vals))
    end
    return vals
end

function basefun_vals_safe!(itp::LinearInterpolation,vals,xs::Vx,x) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    basefun_vals!(itp,vals,xs,x)
end


# function getweight(TT::TransitionTensor{N,k,2,T,probT,pdT},i,j) where {T,probT,pdT<:PDGrid{N,k,T2,NTuple{M,xT}}} where {T2,M,xT<:Axis{itpT}} where {itpT<:LinearInterpolation{true}} where {N,k}
#     return (j == 1 || j == length(TT.pdgrid.xs[i])) ? 0.5 : 1.0
# end
# function getarea(TT::TransitionTensor{N,k,2,T,probT,pdT}) where {T,probT,pdT<:PDGrid{N,k,T2,xsT}} where {T2,xsT<:NTuple{M,aT}} where {M,aT<: Axis{itpT}} where {itpT<:LinearInterpolation{true}} where {N,k}
#     return prod(x.itp.Δ for x in TT.pdgrid.xs) 
# end