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
    # δ³ = δ²*δ

    # ## Assuming explicit derivative in the beginning
    # vals[1] = 0.5δ³ - 0.5δ² - δ + one(δ)
    # vals[2] = -δ³ + δ² + δ

    # ## Midpoint approximated with the same end value
    # # vals[1] = δ³ - 1.5δ² - 0.5δ + one(δ)
    # # vals[2] = -1.5δ³ + 2δ² + 0.5δ

    # ## Mean of the 2 approaches
    # # vals[1] = (2.5δ³ - 3.5δ² - 2δ)/3 + one(δ)
    # # vals[2] = (-4δ³ + 5δ² + 2δ)/3
    # vals[3] = 0.5δ³ - 0.5δ²
    # vals[4] = zero(δ)

    # not enforcing unknown boundary conditions
    vals[1] = 0.5δ² - 1.5δ + one(δ)
    vals[2] = 2δ - δ²
    vals[3] = 0.5δ² - 0.5δ
    vals[4] = zero(δ)
    nothing
end

function cubicinterpolation_weights_end!(vals, dx, Δ)
    δ = dx / Δ;
    δ² = δ*δ
    # δ³ = δ²*δ

    # vals[1] = zero(δ)
    # vals[2] = -0.5δ³ + δ² - 0.5δ
    ## Mean of the 2 approaches
    # vals[3] = (4δ³ - 7δ²)/3 + one(δ)
    # vals[4] = (-2.5δ³ + 4δ²)/3 + 0.5δ

    ## Assuming explicit derivative in the beginning
    # vals[3] = δ³ - 2δ² + one(δ)
    # vals[4] = -0.5δ³ + δ² + 0.5δ

    ## Midpoint approximated with the same end value
    # vals[3] = 1.5δ³ - 2.5δ² + one(δ)
    # vals[4] = -δ³ + 1.5δ² + 0.5δ
    # not enforcing unknown boundary conditions
    vals[1] = zero(δ)
    vals[2] = 0.5δ² - 0.5δ
    vals[4] = one(δ) - δ² 
    vals[3] = 0.5δ² + 0.5δ
    nothing
end

function basefun_vals!(vals,itp::CubicInterpolation,xs::Vx,x) where {Vx<:AbstractVector{Tx}} where Tx<:Number
    basefun_vals_safe!(vals, itp, xs, x)
    nothing
end

function basefun_vals_safe!(vals::SparseInterpolationBaseVals{ord,vT},itp::CubicInterpolation,xs::Vx,x; allow_extrapolation = false, kwargs...) where {ord,vT, Vx<:AbstractVector{Tx}} where Tx<:Number
    needsinterpolation, i = find_idx(xs, x, allow_extrapolation = allow_extrapolation; kwargs...)
    for j in eachindex(vals.val)
        vals.val[j] = zero(eltype(vT))
    end

    if needsinterpolation
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
    elseif i==1
        vals.val[1] = one(eltype(vT))
        vals.val[2] = zero(eltype(vT))
        vals.val[3] = zero(eltype(vT))
        vals.val[4] = zero(eltype(vT))
        vals.idxs[1] = 1
        vals.idxs[2] = 2
        vals.idxs[3] = 3
        vals.idxs[4] = 4
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