function QuinticInterpolation()
    QuinticInterpolation(nothing)
end


function Base.getindex(vals::SparseInterpolationBaseVals{5,vT}, idx) where vT
    if idx == vals.idxs[1]
        return vals.val[1]
    elseif idx == vals.idxs[2]
        return vals.val[2]
    elseif idx == vals.idxs[3]
        return vals.val[3]
    elseif idx == vals.idxs[4]
        return vals.val[4]
    elseif idx == vals.idxs[5]
        return vals.val[5]
    elseif idx == vals.idxs[6]
        return vals.val[6]
    end
    return zero(eltype(vT))
end

function quinticinterpolation_weights!(vals, dx, Δ)
    δ = dx / Δ;
    δ² = δ*δ
    δ³ = δ²*δ
    δ⁴ = δ²*δ²
    δ⁵ = δ²*δ³
    _24= one(δ)/typeof(δ)(24);

    vals[1] = (-5δ⁵ + 13δ⁴ -   δ² + 2δ)*_24 - 0.375δ³
    vals[2] = (25δ⁵ - 64δ⁴ + 16δ² - 16δ)*_24 + 1.625δ³
    vals[3] = (-50δ⁵ - 70δ³)*_24 + 5.25δ⁴ - 1.25δ²+ one(δ)
    vals[4] = (50δ⁵ - 124δ⁴ + 16δ² + 16δ)*_24 + 2.75δ³
    vals[5] = (-25δ⁵ + 61δ⁴ - δ² - 2δ)*_24 - 1.375δ³
    vals[6] = (5δ⁵ + 7δ³)*_24 - 0.5δ⁴

    nothing
end
function quinticinterpolation_weights_1!(vals, dx, Δ) 
    δ = dx / Δ;
    δ² = δ*δ
    δ³ = δ²*δ
    δ⁴ = δ²*δ²
    δ⁵ = δ²*δ³
    _6= one(δ)/typeof(δ)(6);
    
    vals[1] = (7δ⁴ - 4δ³)*_6 - 0.5δ⁵ + 0.5δ² - 1.5δ + one(δ)
    vals[2] = 1.5δ⁵ - 3.5δ⁴ + 2δ³ - δ² + 2δ
    vals[3] = 3.5δ⁴ - 1.5δ⁵ - 2δ³ + 0.5δ² - 0.5δ
    vals[4] = 0.5δ⁵ - (7δ⁴ - 4δ³)*_6

    vals[5] = zero(δ)
    vals[6] = zero(δ)
    nothing
end
function quinticinterpolation_weights_2!(vals, dx, Δ) 
    δ = dx / Δ;
    δ² = δ*δ
    δ³ = δ²*δ
    δ⁴ = δ²*δ²
    δ⁵ = δ²*δ³
    _6= one(δ)/typeof(δ)(6);
    
    vals[1] = (1.25δ⁵ - 2δ)*_6 - 0.5δ⁴ + 0.125δ³ + 0.5δ² 
    vals[2] =  2δ⁴ - (5δ⁵+4δ³)*_6 - δ² - 0.5δ + one(δ)
    vals[3] =  1.25δ⁵ - 3δ⁴ + 1.25δ³ + 0.5δ² + δ
    vals[4] = 2δ⁴ - δ³ - (5δ⁵ + δ)*_6
    vals[5] = (1.25δ⁵ + 1.75δ³)*_6 - 0.5δ⁴

    vals[6] = zero(δ)

    nothing
end
function quinticinterpolation_weights_e1!(vals, dx, Δ) 
    δ = dx / Δ;
    δ² = δ*δ
    δ³ = δ²*δ
    δ⁴ = δ²*δ²
    δ⁵ = δ²*δ³
    _6 = one(δ)/typeof(δ)(6)
    vals[1] = zero(δ)
    vals[2] = zero(δ)
    vals[3] = (8δ⁴ + δ)*_6 - 0.5δ⁵ - δ³
    vals[4] = 1.5δ⁵ - 4δ⁴ + 3δ³ + 0.5δ² - δ
    vals[5] = 4δ⁴ - 1.5δ⁵ - 3δ³ - δ² + 0.5δ + one(δ)
    vals[6] = 0.5δ⁵ + (2δ - 8δ⁴)*_6 + δ³ + 0.5δ²

    nothing
end
function quinticinterpolation_weights_e2!(vals, dx, Δ) 
    δ = dx / Δ;
    δ² = δ*δ
    δ³ = δ²*δ
    δ⁴ = δ²*δ²
    δ⁵ = δ²*δ³
    _24= one(δ)/typeof(δ)(24);
    
    vals[1] = zero(δ)
    vals[2] = (13δ⁴ - 5δ⁵ - δ² + 2δ)*_24 - 0.375δ³
    vals[3] = (20δ⁵ - 52δ⁴ + 32δ³ + 16δ² - 16δ)*_24
    vals[4] = 3.25δ⁴ - 1.25δ⁵ - 1.75δ³ - 1.25δ² + one(δ)
    vals[5] = (20δ⁵ - 52δ⁴ + 16δ² +16δ)*_24 + δ³
    vals[6] = (13δ⁴ - 5δ⁵ - 5δ³ - δ² - 2δ)*_24

    nothing
end

function basefun_vals!(vals,itp::QuinticInterpolation,xs::Vx,x) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    basefun_vals_safe!(vals, itp, xs, x)
    nothing
end

function basefun_vals_safe!(vals::SparseInterpolationBaseVals{ord,vT},itp::QuinticInterpolation,xs::Vx,x; kwargs...) where {ord,vT, Vx<:AbstractVector{Tx}} where Tx<:Number
    do_interpolation, zero_extrapolation, i = find_idx(xs, x; kwargs...)
    # for j in eachindex(vals.val)
    #     vals.val[j] = zero(eltype(vT))
    # end

    if do_interpolation
        if i == 1
            for j in 1:6
                vals.idxs[j] = j
            end
            quinticinterpolation_weights_1!(vals.val, x - xs[i], itp.Δ)
        elseif i == 2
            for j in 1:6
                vals.idxs[j] = j
            end
            quinticinterpolation_weights_2!(vals.val, x - xs[i], itp.Δ)
        elseif i == length(xs) - 1
            vals.idxs[1] = i-4
            vals.idxs[2] = i-3
            vals.idxs[3] = i-2
            vals.idxs[4] = i-1
            vals.idxs[5] = i
            vals.idxs[6] = i+1
            quinticinterpolation_weights_e1!(vals.val, x - xs[i], itp.Δ)
        elseif i == length(xs) - 2
            vals.idxs[1] = i-3
            vals.idxs[2] = i-2
            vals.idxs[3] = i-1
            vals.idxs[4] = i
            vals.idxs[5] = i+1
            vals.idxs[6] = i+2
            quinticinterpolation_weights_e2!(vals.val, x - xs[i], itp.Δ)
        else
            vals.idxs[1] = i-2
            vals.idxs[2] = i-1
            vals.idxs[3] = i
            vals.idxs[4] = i+1
            vals.idxs[5] = i+2
            vals.idxs[6] = i+3
            quinticinterpolation_weights!(vals.val, x - xs[i], itp.Δ)
        end
    elseif zero_extrapolation
        vals.val .= zero(eltype(vT))
    elseif i<7
        @inbounds for j in 1:6
            vals.val[j] = (j == i ? one(eltype(vT)) : zero(eltype(vT)))
            vals.idxs[j] = j
        end
    else
        for j in 1:5
            vals.val[j] = zero(eltype(vT))
        end
        vals.val[6] = one(eltype(vT))
        vals.idxs[1] = i-5
        vals.idxs[2] = i-4
        vals.idxs[3] = i-3
        vals.idxs[4] = i-2
        vals.idxs[5] = i-1
        vals.idxs[6] = i
    end
        
    nothing
    # return vals
end