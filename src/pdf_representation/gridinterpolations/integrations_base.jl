function integrate(f::InterpolatedFunction)
    _integrate(f.p, f.axes...)
end
function _integrate(p::AbstractVector{<:Number},axis::GridAxis)
    sum(axis.wts[i]*_p for (i,_p) in enumerate(p))
end
function _integrate(p::AbstractArray{<:Number,N}, axes::Vararg{Any,N}) where {N}
    sp = size(p,N);
    sum(last(axes).wts[i]*_integrate(view(p,(Colon() for _ in 1:N-1)...,i),axes[1:N-1]) for i in 1:sp)
end

# function expected_value(f::Function, p::InterpolatedFunction)
#     for idx in p.idx_it
#         p.p[idx...]
#     end
# end