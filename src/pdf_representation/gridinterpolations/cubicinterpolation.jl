function CubicInterpolation()
    CubicInterpolation(nothing)
end


function Base.getindex(vals::SparseInterpolationBaseVals{3,vT}, idx) where vT
    if idx == vals.idxs[1]
        return vals.val[1]
    elseif idx == vals.idxs[2]
        return vals.val[2]
    elseif idx == vals.idxs[3]
        return vals.val[3]
    elseif idx == vals.idxs[4]
        return vals.val[4]
    end
    return zero(eltype(vT))
end

function cubicinterpolation_weights!(vals, dx, Δ)
    δ = dx / Δ;
    δ² = δ*δ
    δ³ = δ²*δ

    vals[1] = -0.5δ³ + δ² - 0.5δ
    vals[2] = 1.5δ³ - 2.5δ² + one(δ)
    vals[3] = -1.5δ³ + 2δ² + 0.5δ
    vals[4] = 0.5δ³ - 0.5δ²
    nothing
end
function cubicinterpolation_weights_1!(vals, dx, Δ)
    δ = dx / Δ;
    δ² = δ*δ
    # δ³ = δ²*

    # not enforcing unknown boundary conditions
    # Quadratic
    vals[1] = 0.5δ² - 1.5δ + one(δ)
    vals[2] = 2δ - δ²
    vals[3] = 0.5δ² - 0.5δ
    vals[4] = zero(δ)
    # Cubic (worse)
    # vals[1] = 0.5δ³-0.5δ² - δ + one(δ)
    # vals[2] = δ + δ² - δ³
    # vals[3] = 0.5δ³ - 0.5δ²
    # vals[4] = zero(δ)
    nothing
end

function cubicinterpolation_weights_end!(vals, dx, Δ)
    δ = dx / Δ;
    δ² = δ*δ
    # δ³ = δ²*δ
    
    # not enforcing unknown boundary conditions
    # Quadratic
    vals[1] = zero(δ)
    vals[2] = 0.5δ² - 0.5δ
    vals[3] = one(δ) - δ² 
    vals[4] = 0.5δ² + 0.5δ
    # Cubic (worse)
    # vals[1] = zero(δ)
    # vals[2] = δ² - 0.5δ³ - 0.5δ
    # vals[3] = one(δ) - 2δ² + δ³
    # vals[4] = δ² + 0.5δ - 0.5δ³
    nothing
end

function basefun_vals!(vals,itp::CubicInterpolation,xs::Vx,x) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    basefun_vals_safe!(vals, itp, xs, x)
    nothing
end

function basefun_vals_safe!(vals::SparseInterpolationBaseVals{ord,vT},itp::CubicInterpolation,xs::Vx,x; kwargs...) where {ord,vT, Vx<:AbstractVector{Tx}} where Tx<:Number
    do_interpolation, zero_extrapolation, i = find_idx(xs, x; kwargs...)
    # for j in eachindex(vals.val)
    #     vals.val[j] = zero(eltype(vT))
    # end
    
    if do_interpolation
        if i == 1
            vals.idxs[1] = i
            vals.idxs[2] = i+1
            vals.idxs[3] = i+2
            vals.idxs[4] = i+3
            cubicinterpolation_weights_1!(vals.val, x - xs[i], itp.Δ)
        elseif i == length(xs) - 1
            vals.idxs[1] = i-2
            vals.idxs[2] = i-1
            vals.idxs[3] = i
            vals.idxs[4] = i+1
            cubicinterpolation_weights_end!(vals.val, x - xs[i], itp.Δ)
        else
            vals.idxs[1] = i-1
            vals.idxs[2] = i
            vals.idxs[3] = i+1
            vals.idxs[4] = i+2
            cubicinterpolation_weights!(vals.val, x - xs[i], itp.Δ)
        end
    elseif zero_extrapolation
        vals.val .= zero(eltype(vT))
    elseif i<5
        @inbounds for j in 1:4
            vals.val[j] = (j == i ? one(eltype(vT)) : zero(eltype(vT)))
            vals.idxs[j] = j
        end
    else
        vals.val[1] = zero(eltype(vT))
        vals.val[2] = zero(eltype(vT))
        vals.val[3] = zero(eltype(vT))
        vals.val[4] = one(eltype(vT))
        vals.idxs[1] = i-3
        vals.idxs[2] = i-2
        vals.idxs[3] = i-1
        vals.idxs[4] = i
    end
        
    nothing
    # return vals
end