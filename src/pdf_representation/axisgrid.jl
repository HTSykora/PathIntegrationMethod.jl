Base.getindex(a::AxisGrid,idx...) = a.xs[idx...]
Base.size(a::AxisGrid) = size(a.xs)
_eachindex(axis::AxisGrid) = eachindex(axis)
_eachindex(V::AbstractVector) = eachindex(V)
_gettempvals(axis::AxisGrid) = axis.temp


_eachindex(axis::AxisGrid{itpT}) where itpT<:SparseInterpolationType = _eachindex(axis.temp)
_gettempvals(axis::AxisGrid{itpT}) where itpT<:SparseInterpolationType = axis.temp.val
"""
    LinRange_fromaxis(axis, n)

`n` equidistant points from the first to the last node of `axis`, e.g. to evaluate an [`InterpolatedFunction`](@ref) for plotting.
"""
LinRange_fromaxis(a::AxisGrid, n) = LinRange(a[1],a[end], n)

function sparseinterpolationdata(start, stop, num, xT, wT, newton_cotes_order, order)
    Δ = wT((stop-start)/(num-1));
    xs = LinRange{xT}(start,stop,num)
    wts = NewtonCotesWeights(newton_cotes_order,num,Δ) |> collect
    temp = SparseInterpolationBaseVals(xT,num, order)
    Δ, xs, wts, temp
end

"""
    LinearAxis(start, stop, num; newton_cotes_order = 1, xT = Float64, wT = Float64)

[`AxisGrid`](@ref) of `num` equidistant nodes on `[start, stop]` with piecewise linear interpolation (2 nodes per evaluation).
The quadrature weights are composite Newton–Cotes weights of order `newton_cotes_order` (1: trapezoidal rule, 2: Simpson's rule, 3: Simpson's 3/8 rule;
with end corrections if `num - 1` is not divisible by the order). `xT` and `wT` are the number types of the nodes and of the weights.
"""
function LinearAxis(start,stop,num::Int; xT = Float64, wT = Float64, newton_cotes_order = 1, kwargs...)
    order = 1;
    Δ, xs, wts, temp = sparseinterpolationdata(start, stop, num, xT, wT, newton_cotes_order, order)
    itp = LinearInterpolation(wT(Δ));
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

"""
    CubicAxis(start, stop, num; newton_cotes_order = 2, xT = Float64, wT = Float64)

[`AxisGrid`](@ref) of `num` equidistant nodes on `[start, stop]` with piecewise cubic interpolation: 4 nodes per evaluation
(cubic convolution, i.e. Catmull–Rom splines, and quadratic interpolation in the first and last interval).
For the quadrature weights see [`LinearAxis`](@ref) (default: Simpson's rule).
"""
function CubicAxis(start,stop,num::Int; xT = Float64, wT = Float64, newton_cotes_order = 2, kwargs...)
    order = 3;
    Δ, xs, wts, temp = sparseinterpolationdata(start, stop, num, xT, wT, newton_cotes_order, order)
    itp = CubicInterpolation(wT(Δ));
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

"""
    QuinticAxis(start, stop, num; newton_cotes_order = 2, xT = Float64, wT = Float64)

[`AxisGrid`](@ref) of `num` equidistant nodes on `[start, stop]` with piecewise quintic interpolation: 6 nodes per evaluation
(with modified stencils in the first and last two intervals). For the quadrature weights see [`LinearAxis`](@ref) (default: Simpson's rule).
Recommended for systems with `d > 1`.
"""
function QuinticAxis(start,stop,num::Int; xT = Float64, wT = Float64, newton_cotes_order = 2, kwargs...)
    order = 5;
    Δ, xs, wts, temp = sparseinterpolationdata(start, stop, num, xT, wT, newton_cotes_order, order)
    itp = QuinticInterpolation(wT(Δ));
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

"""
    ChebyshevAxis(start, stop, num; xT = Float64, wT = Float64)

[`AxisGrid`](@ref) of `num` Chebyshev nodes on `[start, stop]` (the extrema of the Chebyshev polynomial of degree `num - 1`, including `start` and `stop`),
with barycentric Lagrange interpolation through every node (dense interpolation) and Clenshaw–Curtis quadrature weights.
The interpolation converges spectrally for smooth PDFs. Recommended for scalar (`d = 1`) systems.
"""
function ChebyshevAxis(start,stop,num::Int; xT = Float64, wT = Float64, kwargs...)
    xs = chebygrid(xT, start,stop,num)
    wts = clenshawcurtisweights(wT, start,stop,num)
    itp = ChebyshevInterpolation(num);
    temp = similar(xs)
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

"""
    TrigonometricAxis(start, stop, num; newton_cotes_order = 2, xT = Float64, wT = Float64)

[`AxisGrid`](@ref) of `num` equidistant nodes on `[start, stop]` with trigonometric (Fourier) interpolation through every node (dense interpolation),
for periodic coordinates (e.g. angles). The period is `stop - start + Δ` with the node distance `Δ = (stop - start)/(num - 1)`, i.e. `stop + Δ` is identified with `start`.
For the quadrature weights see [`LinearAxis`](@ref) (default: Simpson's rule).
"""
function TrigonometricAxis(start,stop,num::Int; xT = Float64, wT = Float64, newton_cotes_order = 2, kwargs...)
    xs = LinRange{xT}(start,stop,num)
    Δ = wT(xs[2] - xs[1])
    wts = NewtonCotesWeights(newton_cotes_order,num,Δ) |> collect
    itp = TrigonometricInterpolation(num,xT(start),xT(stop)+Δ)
    temp = similar(xs)
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

"""
    AxisGrid(start, stop, num; interpolation = :quintic, kwargs...)

Grid of `num` nodes on `[start, stop]` along one coordinate, with an interpolation method and quadrature weights.
The grid of a [`PathIntegration`](@ref) and of an [`InterpolatedFunction`](@ref) is the tensor product of `AxisGrid`s (one for each coordinate).

`interpolation` selects the constructor (`kwargs` are passed on):
- `:linear`: [`LinearAxis`](@ref)
- `:cubic`: [`CubicAxis`](@ref)
- `:quintic`: [`QuinticAxis`](@ref)
- `:chebyshev`: [`ChebyshevAxis`](@ref)
- `:trigonometric`: [`TrigonometricAxis`](@ref)

Sparse interpolations (linear, cubic, quintic) use a few neighbouring nodes, so the step matrix is sparse; dense interpolations (Chebyshev, trigonometric)
use every node of the axis. For `d = 1` [`ChebyshevAxis`](@ref) is recommended, for `d > 1` [`QuinticAxis`](@ref) (or [`CubicAxis`](@ref)).

An `AxisGrid` is an `AbstractVector` of its nodes (`axis[i]`, `axis[end]`, `length(axis)`). Its fields are the interpolation `itp`, the nodes `xs`,
the quadrature weights `wts` and the buffer `temp` of the interpolation weights.
"""
function AxisGrid(start,stop,num::Int; interpolation = :quintic, kwargs...)
    @assert interpolation in [:chebyshev, :linear, :cubic, :quintic, :trigonometric]
    if     interpolation == :linear
        LinearAxis(start, stop, num; kwargs...)
    elseif interpolation == :cubic
        CubicAxis(start, stop, num; kwargs...)
    elseif interpolation == :quintic
        QuinticAxis(start, stop, num; kwargs...)
    elseif interpolation == :chebyshev
        ChebyshevAxis(start, stop, num; kwargs...)
    elseif interpolation == :trigonometric
        TrigonometricAxis(start, stop, num; kwargs...)
    end
end

function remake_gridaxis_with_temp_type(T, ga::AxisGrid{itpT,wT,xT,xeT,tmpT}) where {itpT,wT,xT,xeT,tmpT}
    new_temp = similar(ga.temp, T)
    AxisGrid{itpT,wT,xT,xeT,typeof(new_tmp)}(ga.itp, ga.xs, ga.wts, new_temp)
end

function duplicate(a::AxisGrid{itpT,wT,xT,xeT,tmpT}) where {itpT,wT,xT,xeT,tmpT} 
    AxisGrid{itpT,wT,xT,xeT,tmpT}(deepcopy(a.itp), deepcopy(a.xs), deepcopy(a.wts), deepcopy(a.temp))
end

function duplicate(a::AxisGrid{itpT,wT,xT,xeT,tmpT}) where {itpT<:SparseInterpolationType,wT,xT,xeT,tmpT<:SparseInterpolationBaseVals} 
    AxisGrid{itpT,wT,xT,xeT,tmpT}(deepcopy(a.itp), deepcopy(a.xs), deepcopy(a.wts), zero(a.temp))
end