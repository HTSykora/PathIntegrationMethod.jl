
function DriftTerm(f::Function)
    DriftTerm{1,typeof(f)}(f)
end
function DriftTerm(f::TupleVectorUnion)
    d = length(f)
    DriftTerm{d,typeof(f)}(f)
end
function DriftTerm(f::Vararg{Any,d}) where d
    DriftTerm{d,typeof(f)}(f)
end


function (F::DriftTerm{1,fT})(u,p,t) where fT<:Function
    return F.f(u,p,t)
end
function (F::DriftTerm{1,fT})(i::Integer,u,p,t) where fT<:Function
    return F.f(u,p,t)
end

function (F::DriftTerm{d,fT})(i::Integer,u,p,t) where {d,fT<:TupleVectorUnion}
    if i <= d
        return F.f[i](u,p,t)
    else
        error("i > d")
    end
end

function (F::DriftTerm{d,fT})(du,u,p,t) where {d,fT<:TupleVectorUnion}
    for i in 1:d
        du[i] = D(i,u,p,t)
    end
    return du
end
function (F::DriftTerm{d,fT})(u,p,t) where {d,fT<:TupleVectorUnion}
    du = similar(u) # TODO: type safety!
    F(du,u,p,t)
end