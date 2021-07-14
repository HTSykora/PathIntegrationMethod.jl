Base.getindex(a::Axis,idx...) = a.x[idx...]
Base.size(a::Axis) = size(a.x)

function Axis(start,stop,num::Int; lvl = 1, interpolation = :linear)
    if interpolation == :linear
        xs = LinRange(start,stop,num)
        Δ = (stop-start)/(num-1);
        wts = TrapezoidalWeights{eltype(xs),typeof(Δ)}(num,Δ)
        itp = LinearInterpolation(Δ,Q_equidistant = true);
    elseif interpolation == :chebyshev
        xs = chebygrid(start,stop,num)
        wts = clenshawcurtisweights(start,stop,num)
        itp = ChebyshevInterpolation(num);
    else
        error("Interpolation should be `:linear` or `:chebyshev`")
    end
    return Axis{typeof(itp),eltype(xs),eltype(wts),typeof(xs),typeof(wts)}(itp,xs,wts)
end

function create_temp_axis(T, a::Axis{itpT,xeT,ewT,xT,wT}) where {itpT,xeT,ewT,xT,wT}
    n_itp = new_itp(T, a.itp, length(a))
    Axis{typeof(n_itp),xeT,ewT,xT,wT}(n_itp,a.x,a.wts)
end