function DiscreteIntegrator(discreteintegrator, sdestep::AbstractSDEStep, res_prototype, N::Union{NTuple{1,<:Integer},<:Integer,AbstractArray{<:Integer}}, axes::GA; kwargs...) where GA
    DiscreteIntegrator(discreteintegrator,res_prototype, N, axes; kwargs...)
end

function DiscreteIntegrator(discreteintegrator,res_prototype, N::Union{NTuple{1,<:Integer},<:Integer,AbstractArray{<:Integer}}, axes::GA; xT = Float64, wT = Float64, kwargs...) where GA<:AxisGrid
    start = axes[1]
    stop = axes[end]
    _N = N isa Number ? N : first(N)

    x,w = discreteintegrator(xT, wT, start, stop, _N)
    Q_integrate = Ref(true)
    DiscreteIntegrator{1, typeof(x), typeof(w), typeof(res_prototype), typeof(res_prototype),typeof(Q_integrate)}(x,w,zero(res_prototype),zero(res_prototype),Q_integrate)
end

# QuadGKIntegrator() = QuadGKIntegrator(nothing,nothing,nothing)
QuadGKIntegrator(;kwargs...) = QuadGKIntegrator(nothing,nothing,cleanup_quadgk_keywords(kwargs...))
@inline function cleanup_quadgk_keywords(;σ_init = nothing, μ_init = nothing, allow_extrapolation=false,  zero_extrapolation=true, kwargs...)
    kwargs
end
function DiscreteIntegrator(discreteintegrator::QuadGKIntegrator{T1,T2,Tkwarg},res_prototype, N::Union{NTuple{1,<:Integer},<:Integer}, axes::GA; xT = Float64, wT = Float64, kwargs...) where {GA<:AxisGrid, T1, T2, Tkwarg}
    start = axes[1]
    stop = axes[end]
    # if Tkwargs <: Nothing
    #     qgkkwargs = cleanup_quadgk_keywords(;kwargs...)
    # else
    qgkkwargs = (discreteintegrator.kwargs..., cleanup_quadgk_keywords(;kwargs...)...)
    # end
    QuadGKIntegrator([start, stop], zero(res_prototype), qgkkwargs, Ref(true), zero(res_prototype))
end

function (::ClenshawCurtisIntegrator)(xT, wT, start, stop, num) 
    chebygrid(xT, start, stop, num), clenshawcurtisweights(wT, start, stop, num)
end
function (::GaussLegendreIntegrator)(xT, wT, start, stop, num)
    x0,w0 = gausslegendre(num)
    x = rescale_x(x0, xT, start, stop)
    w = rescale_w(w0, wT, start, stop)
    x, w
end
function (::GaussRadauIntegrator)(xT, wT, start, stop, num)
    x0,w0 = gaussradau(num)
    x = rescale_x(x0, xT, start, stop)
    w = rescale_w(w0, wT, start, stop)
    x, w
end
function (::GaussLobattoIntegrator)(xT, wT, start, stop, num)
    x0,w0 = gausslobatto(num)
    x = rescale_x(x0, xT, start, stop)
    w = rescale_w(w0, wT, start, stop)
    x, w
end

function (::TrapezoidalIntegrator)(xT, wT, start, stop, num) 
    x = LinRange{xT}(start,stop,num)
    Δ = wT((stop-start)/(num-1));
    w = collect(TrapezoidalWeights(num,Δ))
    x, w
end
NewtonCotesIntegrator(N) = NewtonCotesIntegrator{N}()
function (::NewtonCotesIntegrator{N})(xT, wT, start, stop, num)  where N
    x = LinRange{xT}(start,stop,num)
    Δ = wT((stop-start)/(num-1));
    w = collect(NewtonCotesWeights(N,num,Δ))
    x, w
end

function (q::DiscreteIntegrator{nT, xT,wT})(f!, res, temp; Q_reinit_res = true) where {nT, xT<:AbstractVector{T1}, wT<:AbstractVector{T2}} where {T1<:Number, T2<:Number}
    # f!(temp, q.x[1])
    # res .= q.w[1] .* temp
    
    # for (i,x) in enumerate(view(q.x,2:length(q.x)))
    #     f!(temp, x)
    #     res .+= q.w[i+1] .* temp
    # end
    if Q_reinit_res
        res .= zero(eltype(res))
    end
    if q.Q_integrate[]
        for (w,x) in zip(q.w,q.x)
            f!(temp, x)
            temp .*= w
            res .+= temp
        end
    end
end
function (q::DiscreteIntegrator)(f!; kwargs...)
    q(f!, q.res; kwargs...)
end
function (q::DiscreteIntegrator)(f!,res; kwargs...)
    q(f!, res, q.temp; kwargs...)
end

function (q::QuadGKIntegrator)(f!,res; Q_reinit_res = true, kwargs...)
    if !Q_reinit_res
        q.res0 .= res
    end
    if q.Q_integrate
        quadgk!(f!, res, q.int_limits...; q.kwargs...)
    else
        res .= zero(eltype(res))
    end
    if !Q_reinit_res
        res .+= res0
    end
    nothing
end
function (q::QuadGKIntegrator)(f!)
    q(f!, q.res)
end

# Integration schemes for ∫f(x)dx with generic f(x):
# Clenshaw-Curtis (chebyshev with endpoints included)
# Gauss-Legendre (no endpoint is included)
# Gauss-Radau (one endpoint is included)
# Gauss-Lobato (both endpoints are included)
# Newton-Cotes: Trapezoidal vs Simpson rule
# Romberg iteration

# Rescaling
function get_limits(di::QuadGKIntegrator)
    di.int_limits
end
function get_limits(di::DiscreteIntegrator{1})
    (di.x[1], di.x[end])
end
function rescale_to_limits!(di::QuadGKIntegrator,start,stop)
    if isapprox(start,stop, atol = 1.5e-8) || start > stop
        di.Q_integrate[] = false
        return nothing
    end
    di.Q_integrate[] = true
    di.int_limits[1] = start
    di.int_limits[2] = stop
    return nothing
end
function rescale_to_limits!(di::DiscreteIntegrator{1},start,stop)
    if isapprox(start,stop, atol = 1.5e-8) || start > stop
        di.Q_integrate[] = false
        return nothing
    end
    di.Q_integrate[] = true
    rescale_xw!(di.x,di.w,start,stop)
    return nothing
end


function rescale_x(x,T,start,stop)
    bma2 = (T(stop)- T(start))/2
    _start = T(start)
    _x = T.(x)
    _x .= (_x .+ one(T)) .* bma2 .+ _start 
end
function rescale_w(w,T,start,stop)
    bma2 = (T(stop)- T(start))/2
    _x = T.(w)
    _x .= _x .* bma2
end

function rescale_xw!(x,w,start,stop)
    scale = (stop- start)/(x[end] - x[1])
    old_start = x[1];
    x .= (x .- old_start) .* scale .+ start 
    w .= w .* scale
end

function rescale_discreteintegrator!(discreteintegrator::DiscreteIntegrator{1}, sdestep::SDEStep{d,d,m}, pdf; int_limit_thickness_multiplier = 6, kwargs...) where {d,m}
    # discreteintegrator.Q_integrate[] = true
    σ = sqrt(_Δt(sdestep)*get_g(sdestep.sde)(d, sdestep.x0,_par(sdestep),_t0(sdestep))^2)
    mn = min(pdf.axes[d][end], max(pdf.axes[d][1],sdestep.x0[d] - int_limit_thickness_multiplier*σ))
    mx = max(pdf.axes[d][1],min(pdf.axes[d][end],sdestep.x0[d] + int_limit_thickness_multiplier*σ))

    rescale_to_limits!(discreteintegrator,mn,mx)
end