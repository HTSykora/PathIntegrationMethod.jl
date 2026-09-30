Base.getindex(p::InterpolatedFunction,idx...) = p.p[idx...]
Base.size(p::InterpolatedFunction) = size(p.p)
Base.length(p::InterpolatedFunction) = length(p.p)

# Assuming interpolation
"""
    InterpolatedFunction([T = Float64,] axes...; f = nothing)

Function on the tensor product grid of the [`AxisGrid`](@ref)s `axes`, represented by its values `p` at the grid nodes (an `Array{T}` with `size(p) == length.(axes)`).
If `f(x_1, …, x_d)` is given, `p` holds its values at the nodes, otherwise `p` is zero.
The response PDF of a [`PathIntegration`](@ref) (`PI.pdf`) is an `InterpolatedFunction`.

    F(x_1, …, x_d; allow_extrapolation = false, zero_extrapolation = true)

evaluates the interpolation at `(x_1, …, x_d)`. Outside of the grid it returns
- 0 (default),
- the value at the nearest end node of the axis if `zero_extrapolation = false`,
- the interpolation of the first or last interval continued if `allow_extrapolation = true` (periodic continuation for a [`TrigonometricAxis`](@ref)).

`F[i_1, …, i_d]` is the value `F.p[i_1, …, i_d]` at a node. The evaluation uses buffers of the axes: do not evaluate `InterpolatedFunction`s that share an axis
from several threads at the same time.

See also [`integrate`](@ref), [`integrate_diff`](@ref), [`recycle_interpolatedfunction!`](@ref), [`each_latticecoordinate`](@ref).

# Example
```julia
F = InterpolatedFunction(ChebyshevAxis(-1., 1., 21), CubicAxis(0., 2., 41); f = (x, y) -> exp(-x^2 - y))
F(0.3, 1.2)
```
"""
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
    InterpolatedFunction{T,length(psize),itp_type,typeof(axes),typeof(p),typeof(idx_it),typeof(val_it)}(axes,p,idx_it,val_it)
end
InterpolatedFunction(T::DataType,axes::aT; kwargs...) where aT<:Tuple = InterpolatedFunction(T,axes...; kwargs...)

InterpolatedFunction(axes::Vararg{Any,N}; kwargs...) where N = InterpolatedFunction(Float64,axes...; kwargs...)

# Computing interpolated ProbabilityDensityFunction
function (f::InterpolatedFunction{T,N,axesT})(x::Vararg{Any,N}; kwargs...) where {T,N,axesT} #axesT <: NTuple{<:GridAxis}
    interpolate(f.p,f.axes,x...; idx_it = f.idx_it, val_it = f.val_it, kwargs...)
end

function get_itp_type(axes)
    mapreduce(|,axes) do axis
        axis.itp isa SparseInterpolationType 
    end ? SparseInterpolationType : DenseInterpolationType
end
get_val_itp_type(axes) = Val{get_itp_type(axes)}()

is_sparse_interpolation(::InterpolatedFunction{T,N,<:SparseInterpolationType}) where {T,N} = true
is_sparse_interpolation(::InterpolatedFunction) = false
"""
    recycle_interpolatedfunction!(F, f)

Overwrite the node values `F.p` of the [`InterpolatedFunction`](@ref) `F` with the values of the function `f(x_1, …, x_d)` at the grid nodes.
"""
function recycle_interpolatedfunction!(itp_f::InterpolatedFunction, f)
    @assert f isa Function "f is not a function!"
    _it = BI_product(itp_f.axes...)
    for (i,x) in enumerate(_it)
        itp_f.p[i] = f(x...)
    end
end

"""
    each_latticecoordinate(F)

Iterator over the grid nodes `(x_1, …, x_d)` of the [`InterpolatedFunction`](@ref) `F`, in the order of the elements of `F.p` (the first coordinate changes the fastest).

# Example
```julia
F.p .= [exp(-x^2 - v^2) for (x, v) in each_latticecoordinate(F)]
```
"""
function each_latticecoordinate(itp::InterpolatedFunction)
    BI_product(itp.axes...)
end