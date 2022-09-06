# Constructor
function TrigonometricInterpolation(N::Integer, a= 0., b = 2π;)
    c = 2π/(b-a)
    TrigonometricInterpolation{isodd(N),typeof(N),typeof(c)}(N,2π/(b-a))
end

##  Functions for the interpolation
## Odd number of points
function basefun_vals!(vals,itp::TrigonometricInterpolation{true},xs::Vx,x) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    # i = 0... N-1
    Nc = itp.N*itp.c*0.5
    c = itp.c*0.5
    @inbounds for j in 1:itp.N
        vals[j] = sin(Nc*(x-xs[j]))/(itp.N*sin(c*(x-xs[j])))
        # vals[j] = sinc(0.5*Nc*(x-xs[j]))/sinc(0.5*c*(x-xs[j]))
    end
    nothing
    # return vals
end
## Even number of points
function basefun_vals!(vals,itp::TrigonometricInterpolation{false},xs::Vx,x) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    # i = 0... N-1
    Nc = itp.N*itp.c*0.5
    c = itp.c*0.5
    @inbounds for j in 1:itp.N
        vals[j] = sin(Nc*(x-xs[j]))/(itp.N*tan(c*(x-xs[j])))
    end
    # vals[1] *= 0.5
    # vals[end] *= 0.5
    nothing
    # return vals
end

function basefun_vals_safe!(vals,itp::TrigonometricInterpolation,xs::Vx,x; kwargs...) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    do_interpolation, zero_extrapolation, i = find_idx(xs, x; kwargs...)
    if do_interpolation #|| zero(x) < (x - xs[end]) < (xs[2] - xs[1])
        basefun_vals!(vals,itp,xs,x)
    elseif zero_extrapolation
        vals .= zero(eltype(vals))
    else
        vals .= zero(eltype(vals))
        vals[i] = one(eltype(vals))
    end
    nothing
    # return vals
end