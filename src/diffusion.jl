
struct DiffusionTerm{N,k,gT}
    g::gT
end

function DiffusionTerm(g::Function)
    DiffusionTerm{1,1,typeof(g)}(g)
end
function DiffusionTerm(g::Vector,k)
    N = length(g)
    DiffusionTerm{N,k,typeof(g)}(g)
end
function DiffusionTerm(g::Vector)
    DiffusionTerm(g,1)
end

function (D::DiffusionTerm{1,1,gT})(u,p,t) where gT<:Function
    return D.g(u,p,t)
end
function (D::DiffusionTerm{1,1,gT})(n::Integer,u,p,t) where gT<:Function
    return D.g(u,p,t)
end

function (D::DiffusionTerm{N,k,gT})(n::Integer,u,p,t) where {N,k,gT<:Function}
    if n<N
        return zero(eltype(u))
    else 
        return D.g(u,p,t)
    end
end

function (D::DiffusionTerm{N,k,gT})(n::Integer,u,p,t) where {N,k,gT<:Vector}
    if n<k
        return zero(eltype(u))
    end
    if n ∈ k:N
        return D.g[n-k+1](u,p,t)
    end
end
function (D::DiffusionTerm{N,k,gT})(du::T,u::T,p,t) where {N,k,gT<:Vector} where T
    for i in 1:N
        du[i] = D(i,u,p,t)
    end
    return du
end
function (D::DiffusionTerm{N,k,gT})(u,p,t) where {N,k,gT<:Function}
    du = zero(u)
    du[end] = D(i,u,p,t)
    return du
end
function (D::DiffusionTerm{N,k,gT})(u,p,t) where {N,k,gT<:Vector}
    du = zero(u)
    for i in k:N
        du[i] = D(i,u,p,t)
    end
    return du
end
