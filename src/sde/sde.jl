#############################################
# Some utils
_par(sde::AbstractSDE) = sde.par
_par(sde::SDE_VI_Oscillator1D) = sde.osc1D.par
TupleVectorUnion = Union{Tuple,AbstractVector}

# SDE
SDE() = SDE{0,0,Nothing,Nothing,Nothing}(nothing,nothing,nothing)
SDE(d,k) = SDE{d,k,Nothing,Nothing,Nothing}(nothing,nothing,nothing)
SDE(d) = SDE(d,d)
function SDE(f::fT,g::gT; par=nothing) where {fT<:Function,gT<:Function}
    SDE{1,1,DriftTerm{1,1,fT},DiffusionTerm{1,1,gT},typeof(par)}(DriftTerm{1,1,fT}(f),DiffusionTerm{1,1,gT}(g), par)
end

function SDE(f::fT,g::gT; par=nothing) where {fT<:TupleVectorUnion,gT<:TupleVectorUnion} 
    d = length(f); kN = length(g);
    @assert d >= kN "Diffusion term `g` has higher dimensionality then drift term `f`"
    SDE{d,d-kN+1,DriftTerm{d,d-kN+1,fT},DiffusionTerm{d,d-kN+1,gT},typeof(par)}(DriftTerm{d,d-kN+1,typeof(f)}(f),DiffusionTerm{d,d-kN+1,typeof(g)}(g), par)
end

function SDE(N::Integer,f::fT,g::gT; par=nothing) where {fT<:TupleVectorUnion,gT<:Function}
    SDE{N,N,DriftTerm{N,N,fT},DiffusionTerm{N,N,gT},typeof(par)}(DriftTerm{N,N,typeof(f)}(f),DiffusionTerm{N,N,typeof(g)}(g), par)
end

#############################################
# SDE_Oscillator1D
function SDE_Oscillator1D(f::fT,g::gT; par=nothing) where {fT<:Function,gT<:Function}
    SDE_Oscillator1D{DriftTerm{1,fT},DiffusionTerm{1,1,1,1,gT},typeof(par)}(DriftTerm{1,fT}(f),DiffusionTerm{1,1,1,1,gT}(g), par)
end

#############################################
# SDE_VI_Oscillator1D
function SDE_VI_Oscillator1D(f::fT,g::gT, w::Union{Vector{wT},Tuple{wT1,wT2}}; par=nothing) where {fT<:Function,gT<:Function, wT<:Wall, wT1<:Wall, wT2<:Wall}
    sde = SDE_Oscillator1D(f,g, par = par);
    SDE_Oscillator1D{typeof(w), typeof(sde)}(sde, w)
end

SDE_VI_Oscillator1D(f::fT,g::gT, w::Wall; kwargs...) where {fT<:Function,gT<:Function} = SDE_VI_Oscillator1D(f, g, (w,); kwargs...)