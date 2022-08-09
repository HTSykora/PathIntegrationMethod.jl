(w::Wall{<:Function})(v) = w.r(abs(v))
(w::Wall{<:Number})(v) = w.r
Wall(r) = Wall(r,0.);

get_dkm(sde::SDE_VIO) = (2,2,1) # (d,k,m)

# is_v_compatible(v, ID; v_tol = 1.5e-8) = abs(v)<v_tol || !xor(signbit(v), ID == 1)
# is_x_inimpactzone(x, xi, ID; wall_tol = 1.5e-8, kwargs...) = 
# isapprox(x, xi, atol = wall_tol) ? false : xor(ID == 1, x<xi)

osc_f1(u,p,t) = u[2]
function SDE_VIO(f::fT,g::gT, wall::Union{Tuple{wT1},Tuple{wT1,wT2}}, par=nothing) where {fT<:Function,gT<:Function, wT1<:Wall, wT2<:Wall}
    sde = SDE((osc_f1,f), g, par);
    ID = Ref(1)
    SDE_VIO{length(wall),typeof(sde), typeof(wall), typeof(ID)}(sde, wall, ID)
end
SDE_VIO(f::fT,g::gT, w::Wall, par = nothing; kwargs...) where {fT<:Function,gT<:Function} = SDE_VIO(f, g, (w,), par; kwargs...)
