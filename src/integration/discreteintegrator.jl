
function DiscreteIntegrator(discreteintegrator,axes::NTuple{1,GA}, N::NTuple{1,Int},res_prototype; xT = Float64, wT = Float64, kwargs...) where GA<:GridAxis
    start = first(axes)[1]
    stop = first(axes)[end]
    x,w = discreteintegrator(xT, wT, start, stop, first(N))
    DiscreteIntegrator{typeof(discreteintegrator), typeof(x), typeof(w), typeof(res_prototype), typeof(res_prototype)}(x,w,similar(res_prototype),similar(res_prototype))
end

QuadGKIntegrator() = QuadGKIntegrator(nothing,nothing,nothing)
@inline function cleanup_quadgk_keywords(;σ_init = nothing, μ_init = nothing, kwargs...)
    (;kwargs...)
end
function DiscreteIntegrator(discreteintegrator::QuadGKIntegrator, axes::NTuple{1,GA}, N::NTuple{1,Int},res_prototype; xT = Float64, wT = Float64, kwargs...) where GA<:GridAxis
    start = first(axes)[1]
    stop = first(axes)[end]
    QuadGKIntegrator((start, stop),similar(res_prototype),cleanup_quadgk_keywords(;kwargs...))
end

function (::ClenshawCurtisIntegrator)(xT, wT, start, stop, num) 
    chebygrid(xT, start, stop, num), clenshawcurtisweights(wT, num)
end
function (::TrapezoidalIntegrator)(xT, wT, start, stop, num) 
    x = LinRange{xT}(start,stop,num)
    Δ = wT((stop-start)/(num-1));
    w = collect(TrapezoidalWeights{wT}(num,Δ))
    x, w
end


function (q::DiscreteIntegrator{xT,wT})(f!, res, temp) where {xT<:AbstractVector{T1}, wT<:AbstractVector{T2}} where {T1<:Number, T2<:Number}
    f!(temp, q.x[1])
    res .= q.w[1] .* temp
    
    for (i,x) in enumerate(view(q.x,2:length(q.x)))
        f!(temp, x)
        res .= q.w[i+1] .* temp
    end
end
function (q::DiscreteIntegrator)(f!)
    q(f!, q.res)
end
function (q::DiscreteIntegrator)(f!,res)
    q(f!, res, q.temp)
end

function (q::QuadGKIntegrator)(f!,res)
    quadgk!(f!,res,q.int_limits...; q.kwargs...)
end
function (q::QuadGKIntegrator)(f!)
    q(f!,q.res)
end
# Integration schemes for ∫f(x)dx with generic f(x):
# Clenshaw-Curtis (chebyshev with endpoints included)
# Gauss-Legendre (no endpoint is included)
# Gauss-Radau (one endpoint is included)
# Gauss-Lobato (both endpoints are included)
# Romberg iteration
# Newton-Cotes: Trapezoidal vs Simpson rule