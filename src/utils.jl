reduce_tempprod(head::Tuple) = get_tempval(head...)
function reduce_tempprod(head, tail::Vararg{Any,N}) where N
    val = get_tempval(head...)
    val * reduce_tempprod(tail...)
end

BI_product = Base.Iterators.product