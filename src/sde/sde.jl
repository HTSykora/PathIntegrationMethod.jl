#############################################
# Some utils
_par(sde::AbstractSDE) = sde.par
_par(sde::SDE_VIO) = _par(sde.sde)
get_f(sde::AbstractSDE) = sde.f
get_f(sde::SDE_VIO) = get_f(sde.sde)
get_g(sde::AbstractSDE) = sde.g
get_g(sde::SDE_VIO) = get_g(sde.sde)

get_dkm(::AbstractSDE{d,k,m}) where {d,k,m} = (d,k,m)
Q_compatible(::AbstractSDE, ::AbstractSDE) = false
Q_compatible(::AbstractSDE{d,k,m}, ::AbstractSDE{d,k,m}) where {d,k,m}= true
# SDE
SDE() = SDE{0,0,Nothing,Nothing,Nothing}(nothing,nothing,nothing)
SDE(d::Integer,k::Integer) = SDE{d,k,Nothing,Nothing,Nothing}(nothing,nothing,nothing)
SDE(d::Integer) = SDE(d,d)
"""
    SDE(f, g, par = nothing)

Stochastic differential equation of a `d`-dimensional system, where the last coordinates ``x_k, …, x_d`` are driven by independent
Wiener processes ``W_k(t), …, W_d(t)`` (diagonal noise):

    dxᵢ = fᵢ(x, p, t) dt,                        i = 1, …, k-1
    dxᵢ = fᵢ(x, p, t) dt + gᵢ(x, p, t) dWᵢ(t),    i = k, …, d

# Arguments
- `f`: the drift. A function `f(x, p, t)` for a scalar SDE (`d = 1`), or a `Tuple` (or `Vector`) of functions `(f_1, …, f_d)` for `d > 1`,
  where `f_i(x, p, t)` returns the `i`-th component of the drift.
- `g`: the noise intensity of the last coordinate, a function `g(x, p, t)` (`k = d`), or a `Tuple` (or `Vector`) of functions `(g_k, …, g_d)`,
  the noise intensities of the last `d - k + 1` coordinates. The transitional PDF is then integrated over the `d - k + 1` noisy coordinates
  (by default with a [`GaussHermiteIntegrator`](@ref), see [`PathIntegration`](@ref)).
- `par = nothing`: parameters, passed as `p` to `f` and `g`. Use a mutable container (e.g. a `Vector`) to be able to change them later
  with [`recompute_PI!`](@ref) or [`recompute_stepMX!`](@ref). If `par = nothing`, `f` and `g` must not use `p`.

In `f` and `g`, `x` is the state vector (`x[1]`, …, `x[d]`) and `t` is the time.
The drift is traced with Symbolics.jl to compile the backward time step used in the step matrix computation, so `f` has to consist of operations
that accept symbolic arguments (e.g. no branching on the values of `x`, `p` or `t`), unless the back-tracing of [`ExplicitBacktracing`](@ref) is used.

# Example
Duffing oscillator ``dx = v dt``, ``dv = (-2ζv + x - λx³) dt + σ dW(t)`` with `p = [ζ, λ, σ]`:
```julia
f1(x, p, t) = x[2]
f2(x, p, t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x, p, t) = p[3]
sde = SDE((f1, f2), g2, [0.5, 0.25, 1.0])
```
An oscillator driven by white noise and by an Ornstein–Uhlenbeck process ``y`` (noise on ``v`` and ``y``, `k = 2`):
```julia
h1(x, p, t) = x[2]
h2(x, p, t) = -x[1] - 2p[1]*x[2] + x[3]
h3(x, p, t) = -p[2]*x[3]
sde = SDE((h1, h2, h3), ((x, p, t) -> p[3], (x, p, t) -> p[4]), [0.5, 1.0, 0.5, 0.7])
```

See also [`SDE_VIO`](@ref), [`PathIntegration`](@ref).
"""
function SDE(f::fT,g::gT, par=nothing) where {fT<:Function,gT<:Function}
    _f = DriftTerm(f);
    _g = DiffusionTerm(g);
    SDE{1,1,1,typeof(_f),typeof(_g),typeof(par)}(_f,_g, par)
end

# diagonal noise: an independent Wiener process for each of the kN noisy coordinates (m = kN)
function SDE(f::fT,g::gT, par=nothing) where {fT<:TupleVectorUnion,gT<:TupleVectorUnion} 
    d = length(f); kN = length(g); k = d-kN+1;
    @assert d >= kN "Diffusion term `g` has higher dimensionality then drift term `f`"
    SDE{d,k,kN,DriftTerm{d,fT},DiffusionTerm{d,k,kN,kN,gT},typeof(par)}(DriftTerm{d,fT}(f),DiffusionTerm{d,k,kN,kN,gT}(g), par)
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
