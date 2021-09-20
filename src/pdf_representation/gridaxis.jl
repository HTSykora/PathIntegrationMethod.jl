Base.getindex(a::GridAxis,idx...) = a.xs[idx...]
Base.size(a::GridAxis) = size(a.xs)

GridAxis(start,stop,num; kwargs...) = GridAxis(Float64, start,stop,num; kwargs...)
function GridAxis(T, start,stop,num::Int; interpolation = :chebyshev)
    if interpolation == :linear
        xs = LinRange{T}(start,stop,num)
        Δ = (stop-start)/(num-1);
        wts = TrapezoidalWeights{eltype(xs),typeof(Δ)}(num,Δ)
        itp = LinearInterpolation(Δ);
    elseif interpolation == :chebyshev
        xs = chebygrid(T, start,stop,num)
        wts = clenshawcurtisweights(T, start,stop,num)
        itp = ChebyshevInterpolation(num);
    else
        error("Interpolation should be `:linear` or `:chebyshev`.")
    end
    temp = similar(xs)
    return GridAxis{typeof(itp),typeof(wts),typeof(xs),eltype(xs),typeof(temp)}(itp,xs,wts,temp)
end

function remake_gridaxis_with_temp_type(T, ga::GridAxis{itpT,wT,xT,xeT,tmpT}) where {itpT,wT,xT,xeT,tmpT}
    new_temp = similar(ga.temp, T)
    GridAxis{itpT,wT,xT,xeT,typeof(new_tmp)}(ga.itp, ga.xs, ga.wts, new_temp)
end

