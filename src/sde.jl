abstract type AbstractSDE{N,k} end
struct SDE{N,k,fT,gT,pT} <: AbstractSDE{N,k}
    f::fT
    g::gT
    par::pT
end
function SDE(f::fT,g::gT; par=nothing) where {fT<:Function,gT<:Function}
    SDE{1,1,DriftTerm{1,1,fT},DiffusionTerm{1,1,gT},typeof(par)}(DriftTerm{1,1,fT}(f),DiffusionTerm{1,1,gT}(g), par)
end

function SDE(f::fT,g::gT; par=nothing) where {fT<:Vector,gT<:Vector} 
    N = length(f); kN = length(g);
    @assert N >= kN "Diffusion term `g` has higher dimensionality then drift term `f`"
    SDE{N,N-kN+1,DriftTerm{N,N-kN+1,fT},DiffusionTerm{N,N-kN+1,gT},typeof(par)}(DriftTerm{N,N-kN+1,typeof(f)}(f),DiffusionTerm{N,N-kN+1,typeof(g)}(g), par)
end

function SDE(N::Integer,f::fT,g::gT; par=nothing) where {fT<:Vector,gT<:Function}
    SDE{N,N,DriftTerm{N,N,fT},DiffusionTerm{N,N,gT},typeof(par)}(DriftTerm{N,N,typeof(f)}(f),DiffusionTerm{N,N,typeof(g)}(g), par)
end

struct SDE_Oscillator1D{fT, gT, parT} <: AbstractSDE{2,2}
    f::fT
    g::gT
    par::parT
end

function SDE_Oscillator1D(f::fT,g::gT; par=nothing) where {fT<:Function,gT<:Function}
    SDE_Oscillator1D{DriftTerm{1,1,fT},DiffusionTerm{1,1,gT},typeof(par)}(DriftTerm{1,1,fT}(f),DiffusionTerm{1,1,gT}(g), par)
end