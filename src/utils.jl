function (f::fT)(x) where fT <: Scalar_Or_Function{rT} where rT<:Number
    f.r
end
function (f::fT)(x) where fT <: Scalar_Or_Function{rT} where rT<:Function
    f.r(x)
end