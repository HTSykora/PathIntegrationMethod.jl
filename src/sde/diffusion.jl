# g: Rᵈ × [0,T] ↦ Rᵈˣᵐ, gᵢ,ⱼ = 0 for i = 1,...,d-k+1; j = 1 ... m
# and 
# gᵢ,ⱼ ≠ 0 for i = k,...,d; j = 1 ... m
function DiffusionTerm(g::Function)
    DiffusionTerm{1,1,1,1,typeof(g)}(g)
end
function DiffusionTerm(m,g::Function)
    DiffusionTerm{1,1,1,m,typeof(g)}(g)
end
function DiffusionTerm(d,k,m, g::TupleVectorUnion)
    dk = length(g)
    @assert d - k + 1 == dk "Wrong diffusion length: $(d - k + 1) ≠ $(dk)!"
    DiffusionTerm{d,k,dk,m,typeof(g)}(g)
end
function DiffusionTerm(d,k,m, g::Vararg{Any,dk}) where dk
    DiffusionTerm(d,k,m,g)
end

function (D::DiffusionTerm{1,1,1,m,gT})(u,p,t) where {m,gT<:Function}
    return D.g(u,p,t)
end
function (D::DiffusionTerm{1,1,1,m,gT})(i::Integer,u,p,t) where {m,gT<:TupleVectorUnion}
    return D.g(u,p,t)
end

function (D::DiffusionTerm{d,k,dk,m,gT})(i::Integer,u,p,t) where {d,k,dk,m,gT<:Function}
    if i<d
        return zeros(m)
    else 
        return D.g(u,p,t)
    end
end

function (D::DiffusionTerm{d,k,dk,m,gT})(i::Integer,u,p,t) where {d,k,dk,m,gT<:TupleVectorUnion}
    if i<k
        return zeros(m) # TODO: not the same type as the output of g!
    end
    if i ∈ k:d
        return D.g[i-k+1](u,p,t)
    end
end
# function (D::DiffusionTerm{d,k,dk,1,gT})(du,u::T,p,t) where {d,k,dk,gT<:TupleVectorUnion} where T
#     for i in k:d
#         du[i] = D(i,u,p,t)
#     end
#     return du
# end
# function (D::DiffusionTerm{d,k,dk,m,gT})(u,p,t) where {d,k,dk,m,gT<:Function}
#     du = [D(i,u,p,t) for i in 1:d]
#     return du
# end
