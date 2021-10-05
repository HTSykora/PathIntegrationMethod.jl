Base.getindex(a::GridAxis,idx...) = a.xs[idx...]
Base.size(a::GridAxis) = size(a.xs)
_eachindex(axis::GridAxis) = eachindex(axis)
_eachindex(V::AbstractVector) = eachindex(V)
_gettempvals(axis::GridAxis) = axis.temp


_eachindex(axis::GridAxis{itpT}) where itpT<:SparseInterpolationType = _eachindex(axis.temp)
_gettempvals(axis::GridAxis{itpT}) where itpT<:SparseInterpolationType = axis.temp.val

function GridAxis(start,stop,num::Int; xT = Float64, wT = Float64, interpolation = :chebyshev, newton_cotes_order = 2; kwargs...)
    @assert interpolation in [:chebyshev, :linear, :cubic]
    if interpolation == :chebyshev
        xs = chebygrid(xT, start,stop,num)
        wts = clenshawcurtisweights(wT, start,stop,num)
        itp = ChebyshevInterpolation(num);
        temp = similar(xs)
    else 
        Δ = wT((stop-start)/(num-1));
        xs = LinRange{xT}(start,stop,num)
        if interpolation == :linear
            itp = LinearInterpolation(wT(Δ));
            order = 1
        elseif interpolation == :cubic
            itp = CubicInterpolation(wT(Δ));
            order = 3
        end
        wts = NewtonCotesWeights(newton_cotes_order,num,Δ) |> collect
        temp = SparseInterpolationBaseVals(xT,num, order)
    end
    return GridAxis{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

function remake_gridaxis_with_temp_type(T, ga::GridAxis{itpT,wT,xT,xeT,tmpT}) where {itpT,wT,xT,xeT,tmpT}
    new_temp = similar(ga.temp, T)
    GridAxis{itpT,wT,xT,xeT,typeof(new_tmp)}(ga.itp, ga.xs, ga.wts, new_temp)
end

function duplicate(a::GridAxis{itpT,wT,xT,xeT,tmpT}) where {itpT,wT,xT,xeT,tmpT} 
    GridAxis{itpT,wT,xT,xeT,tmpT}(deepcopy(a.itp), deepcopy(a.xs), deepcopy(a.wts), deepcopy(a.temp))
end

function duplicate(a::GridAxis{itpT,wT,xT,xeT,tmpT}) where {itpT<:SparseInterpolationType,wT,xT,xeT,tmpT<:SparseInterpolationBaseVals} 
    GridAxis{itpT,wT,xT,xeT,tmpT}(deepcopy(a.itp), deepcopy(a.xs), deepcopy(a.wts), zero(a.temp))
end