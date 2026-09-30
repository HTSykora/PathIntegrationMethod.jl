
"""
    DriftTerm(f)
    DriftTerm(f_1, …, f_d)

Drift ``f(x, p, t)`` of an [`SDE`](@ref), constructed by `SDE` from a function (`d = 1`) or from the functions of the components
(given as a `Tuple`, a `Vector` or as separate arguments). For a drift `F`:
- `F(i, x, p, t)`: the `i`-th component
- `F(dx, x, p, t)`: all components, written into `dx`
- `F(x, p, t)`: all components (for `d = 1`: the scalar drift)
"""
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
        du[i] = F(i,u,p,t)
    end
    return du
end
function (F::DriftTerm{d,fT})(u,p,t) where {d,fT<:TupleVectorUnion}
    du = similar(u) # TODO: type safety!
    F(du,u,p,t)
end