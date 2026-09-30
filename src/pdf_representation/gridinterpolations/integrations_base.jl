# Product of the quadrature weights of a grid node, w₁⋯w_d (`prod(w)` inside `sum` is not inferred: nested reductions)
node_weight(w::Tuple) = *(w...)

"""
    integrate(F)
    integrate(f, F)

``∫F(x) dx`` over the grid of the [`InterpolatedFunction`](@ref) `F`, computed with the quadrature weights of its axes.
With a function `f(x_1, …, x_d)`: ``∫f(x) F(x) dx``, e.g. an expected value if `F` is a PDF.

# Example
```julia
integrate(PI.pdf)                  # ≈ 1
integrate((x, v) -> x^2, PI.pdf)   # E[x²]
```
"""
function integrate(f::InterpolatedFunction)
    _integrate(f.p, f.axes...)
end
# function _integrate(p::AbstractVector{<:Number},axis::AxisGrid)
#     sum(axis.wts[i]*_p for (i,_p) in enumerate(p))
# end
function _integrate(p::AbstractArray{<:Number,N}, axes::Vararg{Any,N}) where {N}
    # xs = Iterators.product((ax.xs for ax in p.axes)...)
    ws = Iterators.product(map(ax -> ax.wts, axes)...)
    sum(node_weight(w)*p[i] for (i,w) in enumerate(ws))
end
function integrate(f::Function, p::InterpolatedFunction)
    xs = Iterators.product(map(ax -> ax.xs, p.axes)...)
    ws = Iterators.product(map(ax -> ax.wts, p.axes)...)
    sum(f(x...)*node_weight(w)*p.p[i] for (i,(x,w)) in enumerate(zip(xs,ws)))
end 

"""
    integrate_diff(F1, F2; f = abs)
    integrate_diff(F1, p2; f = abs)

``∫f(F₁(x) - F₂(x)) dx`` over the grid of the [`InterpolatedFunction`](@ref) `F1`, computed with the quadrature weights of its axes.
`F2` has to be defined on the same grid; alternatively the array `p2` of its node values can be given. The default `f = abs` gives the ``L¹`` distance.
"""
function integrate_diff(f1::fT, f2::fT; kwargs...) where fT<:InterpolatedFunction
    _integrate_diff(f1.p, f2.p, f1.axes...; kwargs...)
end
function integrate_diff(f1::fT, p2::pT; kwargs...) where {fT<:InterpolatedFunction{T,N},pT<:AbstractArray{T,N}} where {T,N}
    _integrate_diff(f1.p, p2, f1.axes...; kwargs...)
end
function _integrate_diff(p1::AbstractArray{<:Number,N}, p2::AbstractArray{<:Number,N}, axes::Vararg{Any,N}; f = abs, kwargs...) where {N}
    ws = Iterators.product(map(ax -> ax.wts, axes)...)
    sum(node_weight(w)*f(p1[i] - p2[i]) for (i,w) in enumerate(ws))
end



# function expected_value(f::Function, p::InterpolatedFunction)
#     for idx in p.idx_it
#         p.p[idx...]
#     end
# end


# Weight functions
# Base.size(tw::TrapezoidalWeights) = (tw.l,)
# function Base.getindex(tw::TrapezoidalWeights{wT}, idx) where wT
#     if idx == 1 || idx == tw.l
#         return wT(0.5)*tw.Δ
#     else
#         return tw.Δ
#     end
# end

TrapezoidalWeights(l,Δ) = NewtonCotesWeights(1,l,Δ)

NewtonCotesWeights(N, l::Integer, Δ::ΔT) where {lT,ΔT} = NewtonCotesWeights{N,ΔT, mod(l-1,N), typeof(l)}(l, Δ)
Base.size(tw::NewtonCotesWeights) = (tw.l,)
function Base.getindex(ncw::NewtonCotesWeights{1,wT},idx) where wT<:Number
    if idx == 1 || idx == ncw.l
        return wT(0.5)*ncw.Δ
    else
        return ncw.Δ
    end
end

function Base.getindex(ncw::NewtonCotesWeights{2,wT,0},idx) where wT<:Number
    if idx == 1 || idx == ncw.l
        return ncw.Δ/wT(3)
    elseif isodd(idx)
        return wT(2)*ncw.Δ/wT(3)
    else
        return wT(4)*ncw.Δ/wT(3)
    end
end
function Base.getindex(ncw::NewtonCotesWeights{2,wT,1},idx) where wT<:Number
    if idx == 1
        return ncw.Δ / wT(3)
    elseif idx < ncw.l - 3
        if isodd(idx)
            return wT(2) * ncw.Δ/wT(3)
        else
            return wT(4) * ncw.Δ/wT(3)
        end
    elseif idx == ncw.l - 3
        return ncw.Δ * wT(17)/wT(24)
    elseif idx == ncw.l
        return ncw.Δ * wT(0.375)
    else
        return ncw.Δ * wT(1.125)
    end
end

function Base.getindex(ncw::NewtonCotesWeights{3,wT,0},idx) where wT<:Number
    if idx == 1 || idx == ncw.l
        return ncw.Δ * wT(0.375)
    elseif mod(idx - 1, 3) == 0
        return ncw.Δ * wT(0.75)
    else
        return ncw.Δ * wT(1.125)
    end
end
function Base.getindex(ncw::NewtonCotesWeights{3,wT,1},idx) where wT<:Number
    if idx == 1
        return ncw.Δ * wT(0.375)
    elseif idx < ncw.l - 4
        if mod(idx - 1, 3) == 0
            return ncw.Δ * wT(0.75)
        else
            return ncw.Δ * wT(1.125)
        end
    elseif idx == ncw.l - 4
        return ncw.Δ * wT(247)/wT(360)
    elseif idx == ncw.l - 2
        return ncw.Δ * wT(8)/wT(15)
    elseif idx == ncw.l
        return ncw.Δ * wT(14)/wT(45)
    else
        return ncw.Δ * wT(64)/wT(45)
    end
end
function Base.getindex(ncw::NewtonCotesWeights{3,wT,2},idx) where wT<:Number
    if idx == 1
        return ncw.Δ * wT(0.375)
    elseif idx < ncw.l - 5
        if mod(idx - 1, 3) == 0
            return ncw.Δ * wT(0.75)
        else
            return ncw.Δ * wT(1.125)
        end
    elseif idx == ncw.l - 5
        return ncw.Δ * wT(203)/wT(288) # 3/8 + 95/288
    elseif idx == ncw.l - 4 || idx == ncw.l - 1
        return ncw.Δ * wT(125)/wT(96)
    elseif idx == ncw.l
        return ncw.Δ * wT(95)/wT(288)
    else
        return ncw.Δ * wT(125)/wT(144)
    end
end