struct ProbabilityDensityFunction{T,N,axesT,pT} <: AbstractArray{T,N}
    axes::axesT
    p::pT
end

# Assuming interpolation
function ProbabilityDensityFunction(T::DataType, axes::Vararg{Any,N}) where N
    psize = length.(axes)
    p = zeros(T,psize...)
    ProbabilityDensityFunction{eltype(T),length(psize),typeof(axes),typeof(p)}(axes,p)
end
ProbabilityDensityFunction(T,axes::aT) where aT<:Tuple = ProbabilityDensityFunction(T,axes...)

ProbabilityDensityFunction(axes::Vararg{Any,N}) = ProbabilityDensityFunction(Float64,axes...)

# Computing interpolated ProbabilityDensityFunction
function (f::ProbabilityDensityFunction{T,N,axesT})(x::Vararg{Any,N}...) where {T,N,axesT} #axesT <: NTuple{<:GridAxis}
    interpolate(f.p,f.axes,x...)
end