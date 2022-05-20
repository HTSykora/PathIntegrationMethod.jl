Base.getindex(a::AxisGrid,idx...) = a.xs[idx...]
Base.size(a::AxisGrid) = size(a.xs)
_eachindex(axis::AxisGrid) = eachindex(axis)
_eachindex(V::AbstractVector) = eachindex(V)
_gettempvals(axis::AxisGrid) = axis.temp


_eachindex(axis::AxisGrid{itpT}) where itpT<:SparseInterpolationType = _eachindex(axis.temp)
_gettempvals(axis::AxisGrid{itpT}) where itpT<:SparseInterpolationType = axis.temp.val
LinRange_fromaxis(a::AxisGrid, n) = LinRange(a[1],a[end], n)

function sparseinterpolationdata(start, stop, num, xT, wT, newton_cotes_order, order)
    Δ = wT((stop-start)/(num-1));
    xs = LinRange{xT}(start,stop,num)
    wts = NewtonCotesWeights(newton_cotes_order,num,Δ) |> collect
    temp = SparseInterpolationBaseVals(xT,num, order)
    Δ, xs, wts, temp
end

function LinearAxis(start,stop,num::Int; xT = Float64, wT = Float64, newton_cotes_order = 1, kwargs...)
    order = 1;
    Δ, xs, wts, temp = sparseinterpolationdata(start, stop, num, xT, wT, newton_cotes_order, order)
    itp = LinearInterpolation(wT(Δ));
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

function CubicAxis(start,stop,num::Int; xT = Float64, wT = Float64, newton_cotes_order = 2, kwargs...)
    order = 3;
    Δ, xs, wts, temp = sparseinterpolationdata(start, stop, num, xT, wT, newton_cotes_order, order)
    itp = CubicInterpolation(wT(Δ));
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

function QuinticAxis(start,stop,num::Int; xT = Float64, wT = Float64, newton_cotes_order = 2, kwargs...)
    order = 5;
    Δ, xs, wts, temp = sparseinterpolationdata(start, stop, num, xT, wT, newton_cotes_order, order)
    itp = QuinticInterpolation(wT(Δ));
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

function ChebyshevAxis(start,stop,num::Int; xT = Float64, wT = Float64, kwargs...)
    xs = chebygrid(xT, start,stop,num)
    wts = clenshawcurtisweights(wT, start,stop,num)
    itp = ChebyshevInterpolation(num);
    temp = similar(xs)
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

function TrigonometricAxis(start,stop,num::Int; xT = Float64, wT = Float64, newton_cotes_order = 2, kwargs...)
    xs = LinRange{xT}(start,stop,num)
    Δ = wT(xs[2] - xs[1])
    wts = NewtonCotesWeights(newton_cotes_order,num,Δ) |> collect
    itp = TrigonometricInterpolation(num,xT(start),xT(stop)+Δ)
    temp = similar(xs)
    return AxisGrid{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

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