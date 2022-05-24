#############################################
# Some utils
_par(sde::AbstractSDE) = sde.par
_par(sde::SDE_VIO) = _par(sde.sde)

get_dkm(::AbstractSDE{d,k,m}) where {d,k,m} = (d,k,m)
Q_compatible(::AbstractSDE, ::AbstractSDE) = false
Q_compatible(::AbstractSDE{d,k,m}, ::AbstractSDE{d,k,m}) where {d,k,m}= true
# SDE
SDE() = SDE{0,0,Nothing,Nothing,Nothing}(nothing,nothing,nothing)
SDE(d::Integer,k::Integer) = SDE{d,k,Nothing,Nothing,Nothing}(nothing,nothing,nothing)
SDE(d::Integer) = SDE(d,d)
function SDE(f::fT,g::gT, par=nothing) where {fT<:Function,gT<:Function}
    _f = DriftTerm(f);
    _g = DiffusionTerm(g);
    SDE{1,1,1,typeof(_f),typeof(_g),typeof(par)}(_f,_g, par)
end

function SDE(f::fT,g::gT, par=nothing; m = 1) where {fT<:TupleVectorUnion,gT<:TupleVectorUnion} 
    d = length(f); kN = length(g); k = d-kN+1;
    @assert d >= kN "Diffusion term `g` has higher dimensionality then drift term `f`"
    SDE{d,k,m,DriftTerm{d,fT},DiffusionTerm{d,k,kN,m,gT},typeof(par)}(DriftTerm{d,fT}(f),DiffusionTerm{d,k,kN,m,gT}(g), par)
end
function SDE(f::fT,g::gT, par=nothing; kwargs...) where {fT<:TupleVectorUnion,gT<:Function} 
    SDE(f,(g,),par; kwargs...)
end


# function SDE(N::Integer,f::fT,g::gT, par=nothing) where {fT<:TupleVectorUnion,gT<:Function}
#     SDE{N,N,DriftTerm{N,N,fT},DiffusionTerm{N,N,gT},typeof(par)}(DriftTerm{N,N,typeof(f)}(f),DiffusionTerm{N,N,typeof(g)}(g), par)
# end

#############################################
# SDE_Oscillator1D
# function SDE_Oscillator1D(f::fT,g::gT, par=nothing) where {fT<:Function,gT<:Function}
#     SDE_Oscillator1D{DriftTerm{1,fT},DiffusionTerm{1,1,1,1,gT},typeof(par)}(DriftTerm{1,fT}(f),DiffusionTerm{1,1,1,1,gT}(g), par)
# end

#############################################
# SDE_VI_Oscillator1D
