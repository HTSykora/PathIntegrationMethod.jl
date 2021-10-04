Base.getindex(p::InterpolatedFunction,idx...) = p.p[idx...]
Base.size(p::InterpolatedFunction) = size(p.p)

# Assuming interpolation
function InterpolatedFunction(T::DataType, axes::Vararg{Any,N}; f = nothing, kwargs...) where N
    psize = length.(axes)
    p = zeros(T,psize...)
    idx_it = BI_product(_eachindex.(axes)...)
    val_it = BI_product(_gettempvals.(axes)...)
    if f isa Function
        _it = BI_product(axes...)
        for (i,x) in enumerate(_it)
            p[i] = f(x...)
        end
    end

    itp_type = get_itp_type(axes)
    InterpolatedFunction{eltype(T),length(psize),itp_type,typeof(axes),typeof(p),typeof(idx_it),typeof(val_it)}(axes,p,idx_it,val_it)
end
InterpolatedFunction(T::DataType,axes::aT; kwargs...) where aT<:Tuple = InterpolatedFunction(T,axes...; kwargs...)

InterpolatedFunction(axes::Vararg{Any,N}; kwargs...) where N = InterpolatedFunction(Float64,axes...; kwargs...)

# Computing interpolated ProbabilityDensityFunction
function (f::InterpolatedFunction{T,N,axesT})(x::Vararg{Any,N}) where {T,N,axesT} #axesT <: NTuple{<:GridAxis}
    interpolate(f.p,f.axes,x...; idx_it = f.idx_it, val_it = f.val_it)
end

function get_itp_type(axes)
    mapreduce(|,axes) do axis
        axis.itp isa SparseInterpolationType 
    end ? SparseInterpolationType : DenseInterpolationType
end