function (f::Scalar_Or_Function{rT})(x) where rT<:Number
    f.f
end
function (f::Scalar_Or_Function{rT})(x) where rT<:Function
    f.f(x)
end

get_f_type(f::Scalar_Or_Function{rT}) where rT<:Number = rT