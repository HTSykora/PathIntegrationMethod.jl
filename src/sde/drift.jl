
function DriftTerm(f::Function)
    DriftTerm{1,1,typeof(f)}(f)
end
function DriftTerm(f::Union{Vector,Tuple},k)
    N = length(f)
    DriftTerm{N,k,typeof(f)}(f)
end
function DriftTerm(f::Union{Vector,Tuple})
    DriftTerm(f,1)
end

function (F::DriftTerm{1,1,fT})(u,p,t) where fT<:Function
    return F.f(u,p,t)
end
function (F::DriftTerm{1,1,fT})(n::Integer,u,p,t) where fT<:Function
    return F.f(u,p,t)
end
function (F::DriftTerm{N,k,fT})(n::Integer,u,p,t) where {N,k,fT<:Function}
    F.f(u,p,t)
end

function (F::DriftTerm{N,k,fT})(n::Integer,u,p,t) where {N,k,fT<:Union{Vector,Tuple}}
    if n <= N
        return F.f[n](u,p,t)
    else
        error("n > N")
    end
end

function (F::DriftTerm{N,k,fT})(du,u,p,t) where {N,k,fT<:Union{Vector,Tuple}}
    for i in 1:N
        du[i] = D(i,u,p,t)
    end
    return du
end
function (F::DriftTerm{N,k,fT})(u,p,t) where {N,k,fT<:Union{Vector,Tuple}}
    du = similar(u)
    for i in 1:N
        du[i] = D(i,u,p,t)
    end
    return du
end