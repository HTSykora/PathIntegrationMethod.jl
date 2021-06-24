abstract type AbstractInterpolationType{N,d,TF} end

struct ChebyshevInterpolation{N,d} <: AbstractInterpolationType{N,d,false} end
struct EquidistantLinearInterpolation{N,d} <: AbstractInterpolationType{N,d,true} end

