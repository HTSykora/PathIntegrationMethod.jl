Base.getindex(p::ProbabilityDensityFunction,idx...) = p.p[idx...]
Base.size(p::ProbabilityDensityFunction) = size(p.p)

# Assuming interpolation
function ProbabilityDensityFunction(T::DataType, axes::Vararg{Any,N}; f = nothing) where N
    psize = length.(axes)
    p = zeros(T,psize...)
    idx_it = Base.Iterators.product(eachindex.(axes)...)
    if f isa Function
        _it = Base.Iterators.product(axes...)
        for (i,x) in enumerate(_it)
            p[i] = f(x...)
        end
    end
    ProbabilityDensityFunction{eltype(T),length(psize),typeof(axes),typeof(p),typeof(idx_it)}(axes,p,idx_it)
end
ProbabilityDensityFunction(T::DataType,axes::aT; kwargs...) where aT<:Tuple = ProbabilityDensityFunction(T,axes...; kwargs...)

ProbabilityDensityFunction(axes::Vararg{Any,N}; kwargs...) where N = ProbabilityDensityFunction(Float64,axes...; kwargs...)

# Computing interpolated ProbabilityDensityFunction
function (f::ProbabilityDensityFunction{T,N,axesT})(x::Vararg{Any,N}) where {T,N,axesT} #axesT <: NTuple{<:GridAxis}
    interpolate(f.p,f.axes,x...; idx_it = f.idx_it)
end