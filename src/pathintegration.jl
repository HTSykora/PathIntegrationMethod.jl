function PathIntegration(sde::AbstractSDE{d,k,m}, method, ts, axes::Vararg{Any,d}; initialise_pdf = true, f = nothing, pre_compute = false, debug_mode = false, kwargs...) where {d,k,m}
    if method isa DiscreteTimeSteppingMethod
        x0 = zeros(d) # * not type safe for autodiff
        x1 = similar(x0)
        sdestep = SDEStep(sde, method, x0, x1, ts; kwargs)
    end
    if initialise_pdf
        if f isa Nothing
            _f = init_DiagonalNormalPDF(axes...; kwargs...)
        else
            _f = f
        end
    else
        _f = nothing
    end
    pdf = InterpolatedFunction(axes...; f = _f, kwargs...)

    step_idx = [0]

    if pre_compute
        # * for CPU parallelisation extend this
        ikt = IK_temp(pdf.idx_it, # = idx_it
                    Tuple(similar(axis.temp) for axis in axes), # = itpVs
                    similar(pdf.p)) # = itpM
        IK = IntegrationKernel(sdestep, nothing, axes[k:end], ts, pdf, ikt, kwargs)
        if debug_mode
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
    PI.stepMX[PI.step_idx[1]]
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

