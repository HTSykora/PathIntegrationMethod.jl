function PathIntegration(sde::AbstractSDE{d,k,m}, method, ts, axes::Vararg{Any,d};
    discreteintegrator = d==k ? QuadGKIntegrator() : ClenshawCurtisIntegrator(),
    di_N = nothing, di_mul = 10, # discrete integration resolution
    initialise_pdf = true, f = nothing, pre_compute = false, debug_mode = Val{false}, kwargs...) where {d,k,m}
    if method isa DiscreteTimeSteppingMethod
        x0 = zeros(d) # * not type safe for autodiff
        x1 = similar(x0)
        sdestep = SDEStep(sde, method, x0, x1, ts; kwargs...)
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

