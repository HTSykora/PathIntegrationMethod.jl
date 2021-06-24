abstract type AbstractInterpolationType{N,d,TF} end

struct ChebyshevInterpolation{N,d} <: AbstractInterpolationType{N,d,false} end
ChebyshevInterpolation(N,d) = ChebyshevInterpolation{N,d}()
ChebyshevInterpolation(N) = ChebyshevInterpolation(N,1)

struct EquidistantLinearInterpolation{N,d} <: AbstractInterpolationType{N,d,true} end
EquidistantLinearInterpolation(N,d) = EquidistantLinearInterpolation{N,d}()
EquidistantLinearInterpolation(N) = EquidistantLinearInterpolation(N,1)

