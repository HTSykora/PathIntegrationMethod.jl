(w::Wall{<:Function})(v) = w.r(abs(v))
(w::Wall{<:Number})(v) = w.r
"""
    Wall(r, pos = 0.0)

Rigid wall at position `pos` for an [`SDE_VIO`](@ref), with the restitution coefficient `r`: a number, or a function `r(|v|)` of the impact speed.
`w(v)` returns the restitution coefficient of the wall `w` for the impact velocity `v`.
"""
Wall(r) = Wall(r,0.);

get_dkm(sde::SDE_VIO) = (2,2,1) # (d,k,m)

# is_v_compatible(v, ID; v_tol = 1.5e-8) = abs(v)<v_tol || !xor(signbit(v), ID == 1)
# is_x_inimpactzone(x, xi, ID; wall_tol = 1.5e-8, kwargs...) = 
# isapprox(x, xi, atol = wall_tol) ? false : xor(ID == 1, x<xi)

osc_f1(u,p,t) = u[2]
"""
    SDE_VIO(f, g, walls, par = nothing)

Single degree of freedom vibro-impact oscillator

    dx = v dt
    dv = f(u, p, t) dt + g(u, p, t) dW(t),    u = [x, v]

with one or two rigid [`Wall`](@ref)s: at an impact the velocity is reversed and multiplied by the restitution coefficient of the wall, ``v⁺ = -r(|v⁻|) v⁻``.

# Arguments
- `f`: the acceleration `f(u, p, t)`
- `g`: the noise intensity `g(u, p, t)` of the velocity
- `walls`: a `Wall` (a lower wall: the motion is at `x ≥ wall.pos`) or a tuple `(lower_wall, upper_wall)`
- `par = nothing`: parameters, passed as `p` to `f` and `g`

The position axis of the [`PathIntegration`](@ref) should start at the lower wall (and end at the upper wall), e.g.
`PathIntegration(sde, RK4(), 0.01, QuinticAxis(wall.pos, 4., 101), QuinticAxis(-3., 3., 101))`.
"""
function SDE_VIO(f::fT,g::gT, wall::Union{Tuple{wT1},Tuple{wT1,wT2}}, par=nothing) where {fT<:Function,gT<:Function, wT1<:Wall, wT2<:Wall}
    sde = SDE((osc_f1,f), g, par);
    ID = Ref(1)
    SDE_VIO{length(wall),typeof(sde), typeof(wall), typeof(ID)}(sde, wall, ID)
end
SDE_VIO(f::fT,g::gT, w::Wall, par = nothing; kwargs...) where {fT<:Function,gT<:Function} = SDE_VIO(f, g, (w,), par; kwargs...)
