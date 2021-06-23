struct SDE{N,k,fT,gT,pT}
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
    SDE{N,N-kN,DriftTerm{N,N-kN,fT},DiffusionTerm{N,N-kN,gT},typeof(par)}(DriftTerm{N,N-kN,typeof(f)}(f),DiffusionTerm{N,N-kN,typeof(gT)}(g), par)
end

function SDE(N::Integer,f::fT,g::gT; par=nothing) where {fT<:Vector,gT<:Function}
    SDE{N,N-1,DriftTerm{N,N-1,fT},DiffusionTerm{N,N-1,gT},typeof(par)}(DriftTerm{N,N-1,typeof(f)}(f),DiffusionTerm{N,N-1,typeof(g)}(g), par)
end