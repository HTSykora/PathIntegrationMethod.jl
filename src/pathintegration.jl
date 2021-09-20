function PathIntegration(sde::AbstractSDE{d,k,m}, method, ts, axes::Vararg{Any,d}, initialise_pdf = true, f = nothing, pre_compute = false; kwargs...) where {d,k,m}
    if method isa DiscreteTimeSteppingMethod
        sdestep = SDEStep(sde, method, ts; kwargs)
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
    pdf = InterpolatedFunction(axes...; f = f, kwargs...)

    step_idx = [0]

    if pre_compute
        # TODO: for CPU parallelisation extend this
        ikt = IK_temp(pdf.idx_it, # = idx_it
                    [similar(axis.temp) for axis in axes], # = itpVs
                    similar(pdf.p)) # = itpM
        IK = IntegrationKernel(sdestep, nothing, axes[k:end], zeros(Int, d), ts, pdf, ikt)
        step_MX = compute_stepMX(IK; kwargs...)
    else
        step_MX = nothing
        IK = nothing
    end

    PathIntegration(sdestep, pdf, ts, step_MX, step_idx, IK, kwargs)
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

(f::DiagonalNormalPDF)(x) = prod(normal1D_σ2(μ, σ²,x) for (μ, σ²) in zip(f.μ, f.σ²))
