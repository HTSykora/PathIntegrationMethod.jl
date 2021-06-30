abstract type AbstractAxis{T} <:AbstractVector{T} end 

struct Axis{itpT,T,xT} <: AbstractAxis{T}
    itp::itpT
    x::xT
end

Base.getindex(a::Axis,idx...) = a.x[idx...]
Base.size(a::Axis) = size(a.x)

function Axis(start,stop,num::Int; lvl = 1, interpolation = :linear)
    if interpolation == :linear
        xs = LinRange(start,stop,num)
        itp = LinearInterpolation(lvl,(stop-start)/(num-1),Q_equidistant = true);
    elseif interpolation == :chebyshev
        xs = chebygrid(start,stop,num)
        itp = ChebyshevInterpolation(num);
    else
        error("Interpolation should be `:linear` or `:chebyshev`")
    end
    return Axis{typeof(itp),eltype(xs),typeof(xs)}(itp,xs)
end

