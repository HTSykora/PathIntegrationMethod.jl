"""
    PathIntegration(sde, method, ts, axes...; kwargs...)
    
Compute a `PathIntegration` object for computing the response probability density function evolution of the stochastic dynamical system defined by `sde`.

# Arguments
- `sde::AbstractSDE{d,k,m}`: ``d``-dimensional dynamic system which is subjected to an ``m``-dimensional Wiener process acting on the ``k...d`` coordinates.
- `method::`
# Keyword Arguments

# Examples
```julia-repl
julia> bar([1, 2], [1, 2])
1
```
"""
function PathIntegration(sde::AbstractSDE{d,k,m}, method, ts, axes::Vararg{Any,d};
    discreteintegrator = d==k ? QuadGKIntegrator() : ClenshawCurtisIntegrator(),
    di_N = nothing, di_mul = 10, # discrete integration resolution
    initialise_pdf = true, f_init = nothing, pre_compute = false, debug_mode = Val{false}, kwargs...) where {d,k,m}
    if method isa DiscreteTimeSteppingMethod
        x0 = zeros(d) # * not type safe for autodiff
        x1 = similar(x0)
        sdestep = SDEStep(sde, method, x0, x1, ts; kwargs...)
    end
    if initialise_pdf
        if f_init isa Nothing
            _f = init_DiagonalNormalPDF(axes...; kwargs...)
        else
            _f = f_init
        end
    else
        _f = nothing
    end
    pdf = InterpolatedFunction(axes...; f = _f, kwargs...)

    step_idx = [0]

    if pre_compute
        # * for CPU parallelisation extend this
        itpVs = Tuple(zero(axis.temp) for axis in axes); # = itpVs;
        
        ikt = IK_temp(BI_product(_eachindex.(itpVs)...), # = idx_it
                        BI_product(_val.(itpVs)...), # = val_it
                        itpVs, # = itpVs
                        zero(pdf.p)) # = itpM
        if di_N isa Nothing
            if di_mul isa Number
                di_res = Tuple(di_mul*length(ax) for ax in axes[k:end])
            elseif length(di_mul) == k-d+1
                di_res = Tuple(di_mul*length(ax) for ax in axes[k:end])
            else
                error("Incorrect di_mul: should be a `Number` or an `Array`/`Tuple` with length $(d-k+1)")
            end
        else
            di_res = di_N
        end
        di = DiscreteIntegrator(discreteintegrator, pdf.p, di_res, axes[k:end]...; kwargs...)
        IK = IntegrationKernel(sdestep, nothing, di, ts, pdf, ikt, kwargs)
        if debug_mode in (Val{true},true)
            return IK, kwargs
        end
        step_MX = compute_stepMX(IK; kwargs...)
    else
        step_MX = nothing
        IK = nothing
    end
    p_temp = similar(pdf.p);
    PathIntegration(sdestep, pdf, p_temp,ts, step_MX, step_idx, IK, kwargs)
end
_val(vals) = vals
# PathIntegration{dynT, pdT, tsT, tpdMX_type, Tstp_idx, IKT, kwargT}
function advance!(PI::PathIntegration)
    mul!(vec(PI.p_temp), next_stepMX(PI), vec(PI.pdf.p))
    PI.pdf.p .= PI.p_temp ./_integrate(PI.p_temp, PI.pdf.axes...)
    nothing
end

@inline function next_stepMX(PI::PathIntegration{dynT, pdT,tsT}) where {dynT, pdT,tsT<:Number}
    PI.step_MX
end
@inline function next_stepMX(PI::PathIntegration{dynT, pdT,tsT}) where {dynT, pdT,tsT<:AbstractArray}
    PI.step_idx[1] =  mod1(PI.step_idx[1] + 1, length(PI.step_MX))
    PI.step_MX[PI.step_idx[1]]
end

# Computations utilites
function init_DiagonalNormalPDF(axes...; μ_init = nothing, σ_init = nothing, kwargs...)
    if μ_init isa Nothing
        μs = [(axis[end]+axis[1])/2 for axis in axes]
    elseif μ_init isa Number
        μs = [μ_init for _ in axes]
    elseif μ_init isa Vector
        @assert length(μ_init) == length(axes) "Wrong number of initial μ values are given"
        μs = μ_init
    end

    if σ_init isa Nothing
        σ²s = [((axis[end]-axis[1])/12)^2 for axis in axes]
    elseif σ_init isa Number
        σ²s = [σ_init^2 for _ in axes]
    elseif σ_init isa Vector
        @assert length(σ_init) == length(axes) "Wrong number of initial σ values are given"
        σ²s = σ_init.^2
    end

    DiagonalNormalPDF(μs, σ²s)
end

(f::DiagonalNormalPDF)(x...) = prod(normal1D_σ2(μ, σ², _x) for (μ, σ², _x) in zip(f.μ, f.σ², x))

(PI::PathIntegration)(x...) = PI.pdf(x...)

## Recompute functions
function reinit_PI_pdf!(PI,f = nothing)
    if f isa Nothing
        _f = init_DiagonalNormalPDF(PI.pdf.axes...)
    elseif f isa Function
        _f = f
    end
    recycle_interpolatedfunction!(PI.pdf, _f)
end


function recompute_step_MX!(PI::PathIntegration; par = nothing, t = nothing, f = nothing, Q_reinit = false)
    if !(par isa Nothing)
        PI.IK.sdestep.sde.par .= par;
    end

    if t isa AbstractVector
        if length(PI.IK.t) != length(t)
            resize!(PI.IK.t,length(t))
        end
        PI.IK.t .= t
    elseif t isa Number
        resize!(PI.IK.t,2);
        PI.IK.t[1] = zero(eltype(PI.IK.t))
        PI.IK.t[2] = t
    end
    if Q_reinit
        reinit_PI_pdf!(PI)
    elseif f isa Function
        reinit_PI_pdf!(PI,f)
    end

    reinit_stepMX!(PI.step_MX)
    fill_stepMX_ts!(PI.step_MX, PI.IK; PI.IK.kwargs...)
    nothing
end

function reinit_stepMX!(step_MX::AbstractMatrix{T}) where T
    fill!(step_MX,zero(T))
end
function reinit_stepMX!(step_MX::AbstractVector{amT}) where amT<:AbstractMatrix{T} where T
    fill!.(step_MX,zero(T))
end