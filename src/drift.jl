
struct DriftTerm{N,k,fT}
    f::fT
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

function (F::DriftTerm{N,k,fT})(n::Integer,u,p,t) where {N,k,fT<:Vector}
    if n <= N
        return F.f[n](u,p,t)
    else
        error("n > N")
    end
end

function (F::DriftTerm{N,k,fT})(du,u,p,t) where {N,k,fT<:Vector}
    for i in 1:N
        du[i] = D(i,u,p,t)
    end
    return du
end
function (F::DriftTerm{N,k,fT})(u,p,t) where {N,k,fT<:Vector}
    du = similar(u)
    for i in 1:N
        du[i] = D(i,u,p,t)
    end
    return du
end