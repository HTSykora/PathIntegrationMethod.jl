get_N(N::Integer, dim) = Tuple(N for _ in 1:dim)
function get_N(N::NTuple{n,Integer}, dim) where n 
    @assert dim == n "ERROR: incompatible integrator dimensions: length(N) != dim!"
    N
end
# Gauss–Legendre for noise on a single coordinate, Gauss–Hermite (with far fewer nodes) for noise on several coordinates
default_di_N(::AbstractSDE{d,k,m}) where {d,k,m} = d - k + 1 == 1 ? 31 : 7
function defaultdiscreteintegrator(sde::AbstractSDE{d,k,m}, di_N = default_di_N(sde)) where {d,k,m}
    d - k + 1 == 1 ? GaussLegendreIntegrator(di_N) : GaussHermiteIntegrator(di_N, dim = d - k + 1)
end
# (the default of an SDE step depends on the back-tracing: see StrangSplitting in backtracing.jl)
default_di_N(sdestep::AbstractSDEStep) = default_di_N(sdestep.sde)
defaultdiscreteintegrator(sdestep::AbstractSDEStep, di_N = default_di_N(sdestep)) = defaultdiscreteintegrator(sdestep.sde, di_N)

getintegration_dimensions(::AbstractDiscreteIntegratorType{n}) where n = n
function DiscreteIntegrator(discreteintegrator, sdestep::AbstractSDEStep, res_prototype, axes::Vararg{Any,n}; kwargs...) where n
    DiscreteIntegrator(discreteintegrator, res_prototype, axes...; kwargs...)
end

"""
    DiscreteIntegrator(method, res_prototype, axis; xT = Float64, wT = Float64)

The quadrature rule `method` (e.g. [`GaussLegendreIntegrator`](@ref)`(31)`) on the interval `[axis[1], axis[end]]`, for array-valued integrands
(for a [`QuadGKIntegrator`](@ref) an adaptive integrator). For a discrete integrator `q`,

    q(f!, res; Q_reinit_res = true)

computes ``res = Σₖ wₖ f(xₖ)``, where `f!(v, x)` writes every element of ``f(x)`` into `v` (an array like `res_prototype`); `Q_reinit_res = false` adds the integral to `res`.

[`PathIntegration`](@ref) builds the discrete integrator from its `discreteintegrator` keyword argument, and (with `smart_integration = true`) rescales the nodes
for every row of the step matrix to the region where the transitional PDF is not negligible.
"""
function DiscreteIntegrator(discreteintegrator::AbstractDiscreteIntegratorMethod{1},res_prototype, axes::GA; xT = Float64, wT = Float64, kwargs...) where GA<:AxisGrid
    x, w = reference_rule(discreteintegrator, xT, wT, axes)
    DiscreteIntegrator{1}(x, w, zero(res_prototype), zero(res_prototype), Ref(true), copy(x), copy(w), rule_type(discreteintegrator))
end
# Noise on several coordinates: the tensor product of the rules (a rule with dim = 1 is used in every noisy coordinate).
# The nodes are SVectors of the integration variables, the weights are the products of the weights.
function DiscreteIntegrator(discreteintegrator::AbstractDiscreteIntegratorMethod{dim}, res_prototype, axes::Vararg{AxisGrid,n}; xT = Float64, wT = Float64, kwargs...) where {dim,n}
    dim == 1 || dim == n || throw(ArgumentError("The discrete integrator has $dim dimensions, but there are $n noisy coordinates"))
    refs = ntuple(j -> reference_rule(one_dimensional(discreteintegrator, dim == 1 ? 1 : j), xT, wT, axes[j]), Val(n))
    x_ref, w_ref = map(first, refs), map(last, refs)
    l = prod(map(length, x_ref))
    di = DiscreteIntegrator{n}(Vector{SVector{n,xT}}(undef, l), Vector{wT}(undef, l), zero(res_prototype), zero(res_prototype), Ref(true), x_ref, w_ref, rule_type(discreteintegrator))
    fill_tensor_nodes!(di, (j, i) -> x_ref[j][i], (j, i) -> w_ref[j][i])
    di
end
# The nodes and the weights of a rule in one dimension: on [axis[1], axis[end]] for the rules on an interval,
# the standard nodes ξ and the weights ω exp(ξ²) for Gauss–Hermite (∫f(z)dz ≈ √2 s Σ ω exp(ξ²) f(c + √2 s ξ), exact for f = Gaussian N(c, s²) times a polynomial)
function reference_rule(discreteintegrator::AbstractDiscreteIntegratorMethod{1}, xT, wT, axis)
    _x, _w = discreteintegrator(xT, wT, axis[1], axis[end])
    collect(_x), collect(_w) # (LinRange nodes cannot be rescaled in place)
end
function reference_rule(discreteintegrator::GaussHermiteIntegrator{1}, xT, wT, axis)
    ξ, ω = gausshermite(first(discreteintegrator.N))
    xT.(ξ), wT.(ω .* exp.(ξ .^ 2))
end
rule_type(::AbstractDiscreteIntegratorMethod) = IntervalRule()
rule_type(::GaussHermiteIntegrator) = GaussianRule()
# The rule of the j-th dimension
one_dimensional(di::ClenshawCurtisIntegrator, j) = ClenshawCurtisIntegrator(di.N[j])
one_dimensional(di::GaussLegendreIntegrator, j) = GaussLegendreIntegrator(di.N[j])
one_dimensional(di::GaussRadauIntegrator, j) = GaussRadauIntegrator(di.N[j])
one_dimensional(di::GaussLobattoIntegrator, j) = GaussLobattoIntegrator(di.N[j])
one_dimensional(di::TrapezoidalIntegrator, j) = TrapezoidalIntegrator(di.N[j])
one_dimensional(di::NewtonCotesIntegrator{dim,ord}, j) where {dim,ord} = NewtonCotesIntegrator(di.N[j], ord)
one_dimensional(di::GaussHermiteIntegrator, j) = GaussHermiteIntegrator(di.N[j])
# The tensor product nodes SVector(node(1, i₁), …, node(n, iₙ)) and weights weight(1, i₁) ⋯ weight(n, iₙ)
function fill_tensor_nodes!(di::DiscreteIntegrator{n}, node, weight) where n
    for (l, I) in enumerate(CartesianIndices(map(length, di.x_ref)))
        di.x[l] = SVector(ntuple(j -> node(j, I[j]), Val(n)))
        di.w[l] = prod(ntuple(j -> weight(j, I[j]), Val(n)))
    end
    di
end

# QuadGKIntegrator(; rtol, atol, maxevals, order...):
"""
    QuadGKIntegrator(; kwargs...)

Adaptive Gauss–Kronrod quadrature with `QuadGK.quadgk!`; the keyword arguments (e.g. `rtol`, `atol`, `maxevals`, `order`) are passed to `quadgk!`.
It is much slower than the fixed quadrature rules, but useful as a reference to check their accuracy. Use it as the `discreteintegrator` of [`PathIntegration`](@ref).
"""
QuadGKIntegrator(;kwargs...) = QuadGKIntegrator(nothing,nothing,cleanup_quadgk_keywords(;kwargs...),nothing,nothing)
@inline function cleanup_quadgk_keywords(;σ_init = nothing, μ_init = nothing, allow_extrapolation=false,  zero_extrapolation=true, kwargs...)
    kwargs
end
function DiscreteIntegrator(discreteintegrator::QuadGKIntegrator, res_prototype, axes::GA; kwargs...) where {GA<:AxisGrid}
    start = axes[1]
    stop = axes[end]
    QuadGKIntegrator([start, stop], zero(res_prototype), discreteintegrator.kwargs, Ref(true), zero(res_prototype))
end
DiscreteIntegrator(::QuadGKIntegrator, res_prototype, axes::Vararg{AxisGrid,n}; kwargs...) where n =
    throw(ArgumentError("QuadGKIntegrator is only available for noise on a single coordinate"))

function (di::ClenshawCurtisIntegrator{1})(xT, wT, start, stop) 
    chebygrid(xT, start, stop, first(di.N)), clenshawcurtisweights(wT, start, stop, first(di.N))
end
function (di::GaussLegendreIntegrator{1})(xT, wT, start, stop)
    x0,w0 = gausslegendre(first(di.N))
    x = rescale_x(x0, xT, start, stop)
    w = rescale_w(w0, wT, start, stop)
    x, w
end
function (di::GaussRadauIntegrator{1})(xT, wT, start, stop)
    x0,w0 = gaussradau(first(di.N))
    x = rescale_x(x0, xT, start, stop)
    w = rescale_w(w0, wT, start, stop)
    x, w
end
function (di::GaussLobattoIntegrator{1})(xT, wT, start, stop)
    x0,w0 = gausslobatto(first(di.N))
    x = rescale_x(x0, xT, start, stop)
    w = rescale_w(w0, wT, start, stop)
    x, w
end

function (di::TrapezoidalIntegrator{1})(xT, wT, start, stop) 
    num = first(di.N)
    x = LinRange{xT}(start,stop,num)
    Δ = wT((stop-start)/(num-1));
    w = collect(TrapezoidalWeights(num,Δ))
    x, w
end
function (di::NewtonCotesIntegrator{1, ord})(xT, wT, start, stop)  where ord
    num = first(di.N)
    x = LinRange{xT}(start,stop,num)
    Δ = wT((stop-start)/(num-1));
    w = collect(NewtonCotesWeights(ord,num,Δ))
    x, w
end

function (q::DiscreteIntegrator)(f!, res, temp; Q_reinit_res = true, kwargs...)
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
        # f! may only write the elements that change (e.g. the interpolation stencil of the previous node), which is only
        # valid if the last f! call on `temp` was from this loop (e.g. not with several integrators sharing the stencil)
        temp .= zero(eltype(temp))
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
    res .= zero(eltype(res))
    if q.Q_integrate[]
        quadgk!((v, x) -> (fill!(v, zero(eltype(v))); f!(v, x)), res, q.int_limits...; q.kwargs...)
    end
    if !Q_reinit_res
        res .+= q.res0
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
function rescale_to_limits!(di::IntervalDiscreteIntegrator{1},start,stop)
    if isapprox(start,stop, atol = 1.5e-8) || start > stop
        di.Q_integrate[] = false
        return nothing
    end
    di.Q_integrate[] = true
    rescale_xw!(di.x,di.w,di.x_ref,di.w_ref,start,stop)
    return nothing
end
rescale_to_limits!(di::IntervalDiscreteIntegrator{1}, starts::Tuple{Any}, stops::Tuple{Any}) = rescale_to_limits!(di, starts[1], stops[1])
# Several noisy coordinates: the window [starts[j], stops[j]] in each dimension
function rescale_to_limits!(di::IntervalDiscreteIntegrator{n}, starts::NTuple{n,Any}, stops::NTuple{n,Any}) where n
    if any(j -> isapprox(starts[j], stops[j], atol = 1.5e-8) || starts[j] > stops[j], 1:n)
        di.Q_integrate[] = false
        return nothing
    end
    di.Q_integrate[] = true
    scales = ntuple(j -> (stops[j] - starts[j]) / (di.x_ref[j][end] - di.x_ref[j][1]), Val(n))
    fill_tensor_nodes!(di, (j, i) -> (di.x_ref[j][i] - di.x_ref[j][1]) * scales[j] + starts[j], (j, i) -> di.w_ref[j][i] * scales[j])
    return nothing
end
# Gauss–Hermite: the nodes c + √2 s ξ and the weights √2 s ω exp(ξ²) for the Gaussian N(c, s²) in each dimension
function rescale_to_gaussian!(di::GaussianDiscreteIntegrator{1}, cs::Tuple{Any}, ss::Tuple{Any})
    c, s = cs[1], ss[1]
    di.Q_integrate[] = s > zero(s)
    di.x .= c .+ sqrt(2) * s .* di.x_ref
    di.w .= sqrt(2) * s .* di.w_ref
    return nothing
end
function rescale_to_gaussian!(di::GaussianDiscreteIntegrator{n}, cs::NTuple{n,Any}, ss::NTuple{n,Any}) where n
    di.Q_integrate[] = all(s -> s > zero(s), ss)
    fill_tensor_nodes!(di, (j, i) -> cs[j] + sqrt(2) * ss[j] * di.x_ref[j][i], (j, i) -> sqrt(2) * ss[j] * di.w_ref[j][i])
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

# Map the reference nodes spanning [x_ref[1], x_ref[end]] to [start, stop]
# (always from the reference, so the nodes do not depend on the previous rescalings)
function rescale_xw!(x,w,x_ref,w_ref,start,stop)
    scale = (stop- start)/(x_ref[end] - x_ref[1])
    ref_start = x_ref[1];
    x .= (x_ref .- ref_start) .* scale .+ start 
    w .= w_ref .* scale
end

# (pdf::pdfT forces specialisation: InterpolatedFunction <: Function)
function rescale_discreteintegrator!(discreteintegrator::IntervalDiscreteIntegrator{1}, sdestep::SDEStep{d,d,m}, pdf::pdfT; kwargs...) where {d,m,pdfT}
    mn, mx = get_rescale_limits(sdestep, pdf; kwargs...)
    rescale_to_limits!(discreteintegrator, mn, mx)
end
function rescale_discreteintegrator!(discreteintegrator::QuadGKIntegrator, sdestep::SDEStep{d,d,m}, pdf::pdfT; kwargs...) where {d,m,pdfT}
    mn, mx = get_rescale_limits(sdestep, pdf; kwargs...)
    rescale_to_limits!(discreteintegrator, mn, mx)
end
# In general: the window (rules on an interval) or the Gaussian (Gauss–Hermite) of the transitional PDF in each noisy coordinate i = k, …, d
# (diagonal noise), around sdestep.x0[k:d]: the start of the step traced back from the grid point, or the grid point itself (explicit back-tracing)
function rescale_discreteintegrator!(di::DiscreteIntegrator{n}, sdestep::SDEStep{d,k,m}, pdf::pdfT; int_limit_thickness_multiplier = 6, kwargs...) where {n,d,k,m,pdfT}
    cs = ntuple(j -> sdestep.x0[k+j-1], Val(n))
    rescale_rule!(di, di.rule, cs, noise_stds(sdestep, Val(n)), ntuple(j -> pdf.axes[k+j-1], Val(n)), int_limit_thickness_multiplier)
end
rescale_rule!(di, ::GaussianRule, cs, ss, axes, multiplier) = rescale_to_gaussian!(di, cs, ss)
function rescale_rule!(di::DiscreteIntegrator{n}, ::IntervalRule, cs, ss, axes, multiplier) where n
    starts = ntuple(j -> min(axes[j][end], max(axes[j][1], cs[j] - multiplier*ss[j])), Val(n))
    stops = ntuple(j -> max(axes[j][1], min(axes[j][end], cs[j] + multiplier*ss[j])), Val(n))
    rescale_to_limits!(di, starts, stops)
end
# The standard deviations √(Δt gᵢ(x₀)²) of the noisy coordinates i = k, …, d
noise_stds(sdestep::SDEStep{d,k,m}, ::Val{n}) where {d,k,m,n} = ntuple(j -> sqrt(_Δt(sdestep)*get_g(sdestep.sde)(k+j-1, sdestep.x0,_par(sdestep),_t0(sdestep))^2), Val(n))

function get_rescale_limits(sdestep::SDEStep{d,d,m}, pdf::pdfT; int_limit_thickness_multiplier = 6, kwargs...) where {d,m,pdfT}
    σ = sqrt(_Δt(sdestep)*get_g(sdestep.sde)(d, sdestep.x0,_par(sdestep),_t0(sdestep))^2)
    mn = min(pdf.axes[d][end], max(pdf.axes[d][1],sdestep.x0[d] - int_limit_thickness_multiplier*σ))
    mx = max(pdf.axes[d][1],min(pdf.axes[d][end],sdestep.x0[d] + int_limit_thickness_multiplier*σ))
    mn, mx
end