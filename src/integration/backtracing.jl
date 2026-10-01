## Back-tracing of the time steps from the grid points: NewtonBacktracing (the symbolic Newton iteration of sdestep.jl and driftstep.jl),
## ExplicitBacktracing and StrangSplitting (explicit drift steps backward in time)

# The step tracer of an SDEStep
step_tracer(::NewtonBacktracing, sde, method, x0, x1, t0, t1) = PreComputeNewtonStep()(sde, method, x0, x1, t0, t1)
function step_tracer(bt::Union{ExplicitBacktracing,StrangSplitting}, sde::AbstractSDE{d,k,m}, method, x0, x1, t0, t1) where {d,k,m}
    sde isa SDE || throw(ArgumentError("$(nameof(typeof(bt)))() is not available for $(nameof(typeof(sde))), use NewtonBacktracing()"))
    (k == d || bt isa ExplicitBacktracing) || throw(ArgumentError("$(nameof(typeof(bt)))() needs the noise on the last coordinate only"))
    explicit_step_tracer(bt)
end
explicit_step_tracer(::ExplicitBacktracing) = ExplicitStepTracer()
explicit_step_tracer(::StrangSplitting) = TransportStepTracer()

# The back-tracing method of an SDE step
backtracing(::SDEStep{d,k,m,sdeT,methodT,<:AbstractSymbolicNewtonStepTracer}) where {d,k,m,sdeT,methodT} = NewtonBacktracing()
backtracing(::SDEStep{d,k,m,sdeT,methodT,ExplicitStepTracer}) where {d,k,m,sdeT,methodT} = ExplicitBacktracing()
backtracing(::SDEStep{d,k,m,sdeT,methodT,TransportStepTracer}) where {d,k,m,sdeT,methodT} = StrangSplitting()
backtracing(::NonSmoothSDEStep) = NewtonBacktracing()

task_copy(t::Union{ExplicitStepTracer,TransportStepTracer,DiffusionStepTracer}) = t

## Explicit drift steps (also with dual numbers)

# Tag of the dual numbers of the back-tracing
struct BacktracingTag end

# The drift at the state x (an SVector)
drift_svector(F, x::SVector{d}, par, t) where d = SVector(ntuple(i -> F(i, x, par, t), Val(d)))

# One drift step of the explicit method from time t with the step h (h < 0: backward in time)
explicit_driftstep(::Euler, F, x, par, t, h) = x + h * drift_svector(F, x, par, t)
# The stages are kᵢ = h f(x + Σₗ a_{i-1,l} kₗ, t + c_{i-1} h) (k₁ = h f(x, t)), and the step is x + Σᵢ bᵢ kᵢ, summed in the same order as in
# _eval_driftstep! and fill_to_x1!. The loop over the stages is unrolled, as the rows of the Butcher tableau have different types.
@generated function explicit_driftstep(rk::RungeKutta{ord,ButcherTableau{aT,bT,cT,_cT}}, F, x, par, t, h) where {ord,aT,bT,cT,_cT}
    ex = quote
        a, b, c = rk.BT.a, rk.BT.b, rk.BT.c
        k1 = h * drift_svector(F, x, par, t)
        ks1 = (k1,)
    end
    for j in 1:length(aT.parameters)
        ks, ksnew = Symbol(:ks, j), Symbol(:ks, j + 1)
        terms = [:(a[$j][$l].weight * $ks[a[$j][$l].idx]) for l in 1:length(aT.parameters[j].parameters)]
        push!(ex.args, :($ksnew = ($ks..., h * drift_svector(F, +(x, $(terms...)), par, t + c[$j] * h))))
    end
    ks = Symbol(:ks, length(aT.parameters) + 1)
    push!(ex.args, :(return +(x, $([:(b[$l].weight * $ks[b[$l].idx]) for l in 1:length(bT.parameters)]...))))
    ex
end

# The start x₀ of the drift step of `step` that ends at x₁ (an SVector): one drift step backward in time (from t₁ to t₀),
# and the absolute Jacobian determinant |det ∂x₀/∂x₁| (computed with dual numbers)
function backward_driftstep(step::SDEStep{d}, x1::SVector{d,T}) where {d,T}
    seed = SVector(ntuple(i -> ForwardDiff.Dual{BacktracingTag}(x1[i], ForwardDiff.Partials(ntuple(j -> T(i == j), Val(d)))), Val(d)))
    x0 = explicit_driftstep(step.method.drift, get_f(step.sde), seed, _par(step), _t1(step), -_Δt(step))
    J = SMatrix{d,d,T}(ntuple(l -> ForwardDiff.partials(x0[(l - 1) % d + 1], (l - 1) ÷ d + 1), Val(d*d)))
    map(ForwardDiff.value, x0), abs(det(J))
end

# The grid point x₁ with z in the noisy (last) coordinate
with_last(x1, z, ::Val{d}) where d = SVector(ntuple(i -> i < d ? x1[i] : z, Val(d)))
# The grid point x₁ with z (a number or an SVector) in the noisy coordinates k, …, d
with_noisy(x1, z, ::Val{d}, ::Val{k}) where {d,k} = SVector(ntuple(i -> i < k ? x1[i] : z[i-k+1], Val(d)))

## Kernel values: the transitional PDF (with the Jacobian determinant) at the integration variable z; the start of the step is written into sdestep.x0

# ExplicitBacktracing: z is the last coordinate after the drift step, x₀ = Φ⁻¹(x₁ with z in the last coordinate)
function kernel_value!(IK::IntegrationKernel{kd,<:SDEStep{d,d,m,sdeT,methodT,ExplicitStepTracer}}, z) where {kd,d,m,sdeT,methodT}
    step = IK.sdestep
    x0, detJ = backward_driftstep(step, with_last(IK.x1, z, Val(d)))
    step.x0 .= x0
    fx = normal1D_σ2(z, _Δt(step) * get_g(step.sde)(d, step.x0, _par(step), _t0(step))^2, IK.x1[d])
    isapprox(fx, zero(fx), atol = 1e-8) ? zero(fx) : fx * detJ
end
# Noise on the coordinates k, …, d: z (an SVector) are the noisy coordinates after the drift step
function kernel_value!(IK::IntegrationKernel{kd,<:SDEStep{d,k,m,sdeT,methodT,ExplicitStepTracer}}, z) where {kd,d,k,m,sdeT,methodT}
    step = IK.sdestep
    x0, detJ = backward_driftstep(step, with_noisy(IK.x1, z, Val(d), Val(k)))
    step.x0 .= x0
    fx = diagonal_gaussian(step, IK.x1, Tuple(z))
    isapprox(fx, zero(fx), atol = 1e-8) ? zero(fx) : fx * detJ
end
# StrangSplitting, transport: the transitional PDF is a Dirac delta at x₀ = Φ⁻¹(x₁) (integrated with a single node of weight 1)
function kernel_value!(IK::IntegrationKernel{kd,<:SDEStep{d,d,m,sdeT,methodT,TransportStepTracer}}, _) where {kd,d,m,sdeT,methodT}
    x0, detJ = backward_driftstep(IK.sdestep, SVector(ntuple(i -> IK.x1[i], Val(d))))
    IK.sdestep.x0 .= x0
    detJ
end
# StrangSplitting, diffusion: the step starts at x₁ with z in the last coordinate
function kernel_value!(IK::IntegrationKernel{kd,<:SDEStep{d,d,m,sdeT,methodT,DiffusionStepTracer}}, z) where {kd,d,m,sdeT,methodT}
    step = IK.sdestep
    step.x0 .= with_last(IK.x1, z, Val(d))
    fx = normal1D_σ2(z, _Δt(step) * get_g(step.sde)(d, step.x0, _par(step), _t0(step))^2, IK.x1[d])
    isapprox(fx, zero(fx), atol = 1e-8) ? zero(fx) : fx
end

# Integration window: centred at the grid point (sdestep.x0 = x₁ here), without the backward Newton step
function rescale_discreteintegrator!(IK::IntegrationKernel{kd,<:SDEStep{d,k,m,sdeT,methodT,<:Union{ExplicitStepTracer,DiffusionStepTracer}}}; kwargs...) where {kd,d,k,m,sdeT,methodT}
    rescale_discreteintegrator!(IK.discreteintegrator, IK.sdestep, IK.pdf; kwargs...)
end
rescale_discreteintegrator!(::IntegrationKernel{kd,<:SDEStep{d,k,m,sdeT,methodT,TransportStepTracer}}; kwargs...) where {kd,d,k,m,sdeT,methodT} = nothing

## Integration kernels of the step matrix computation

function integration_kernel(sdestep::AbstractSDEStep{d,k,m}, discreteintegrator, pdf, ts, generic_row_kernel, ik_kwargs; kwargs...) where {d,k,m}
    Q_generic = generic_row_kernel || Q_generic_row_kernel(discreteintegrator)
    # only the generic row kernel needs full-size buffers in the discrete integrator
    res_prototype = Q_generic ? pdf.p : similar(pdf.p, ntuple(_ -> 0, d))
    di = DiscreteIntegrator(discreteintegrator, sdestep, res_prototype, pdf.axes[k:end]...; kwargs...)
    IntegrationKernel(sdestep, nothing, di, ts, pdf, kernel_temp(pdf, di, Q_generic), ik_kwargs)
end
kernel_temp(pdf, di, Q_generic) = IK_temp(map(axis -> zero(axis.temp), pdf.axes), zero(pdf.p), Q_generic ? GenericRowKernel() : row_kernel(pdf, length(di.x)))

# StrangSplitting: the kernels of the transport and of the diffusion steps
function integration_kernel(sdestep::SDEStep{d,k,m,sdeT,methodT,TransportStepTracer}, discreteintegrator, pdf, ts, generic_row_kernel, ik_kwargs; kwargs...) where {d,k,m,sdeT,methodT}
    # the tolerances of the step matrix are applied to the product of the matrices
    sub_kwargs = merge(ik_kwargs, (sparse_tol = zero(ik_kwargs.sparse_tol), sparse_rtol = zero(ik_kwargs.sparse_rtol)))
    di = single_node_integrator(generic_row_kernel ? pdf.p : similar(pdf.p, ntuple(_ -> 0, d)))
    transport = IntegrationKernel(sdestep, nothing, di, ts, pdf, kernel_temp(pdf, di, generic_row_kernel), sub_kwargs)
    diffusion = integration_kernel(diffusion_step(sdestep), discreteintegrator, pdf, ts, generic_row_kernel, sub_kwargs; kwargs...)
    StrangKernel(transport, diffusion, sdestep, ts, pdf, ik_kwargs)
end
# A discrete integrator with a single node of weight 1
function single_node_integrator(res_prototype)
    x, w = [0.0], [1.0]
    DiscreteIntegrator{1}(x, w, zero(res_prototype), zero(res_prototype), Ref(true), copy(x), copy(w), IntervalRule())
end
# The diffusion step of the same SDE (with its own states and time interval)
function diffusion_step(s::SDEStep{d,k,m}) where {d,k,m}
    t0, t1 = Ref(s.t0[]), Ref(s.t1[])
    SDEStep{d,k,m,typeof(s.sde),typeof(s.method),DiffusionStepTracer,typeof(s.x0),typeof(s.x1),typeof(t0),Nothing,Nothing,Nothing}(
        s.sde, s.method, similar(s.x0), similar(s.x1), t0, t1, DiffusionStepTracer(), nothing, nothing, nothing)
end

## StrangSplitting step matrices: S = D(tₘ → t₁) T(t₀ → t₁) D(t₀ → tₘ), tₘ = (t₀ + t₁)/2

function fill_stepMX_ts!(stepMX::AbstractVector{<:AbstractMatrix}, K::StrangKernel; rowcomputation = SerialRowComputation(), smart_integration = true, kwargs...)
    # the rows of the operators are computed as the rows of the step matrix: the stored matrices are the transposed operators
    Tᵀ, D₁ᵀ, D₂ᵀ = operator_storage(K.transport), operator_storage(K.diffusion), operator_storage(K.diffusion)
    ws_T = RowWorkspace(K.transport, rowcomputation, Tᵀ)
    ws_D = RowWorkspace(K.diffusion, rowcomputation, D₁ᵀ)
    for jₜ in 1:length(K.t)-1
        t0, t1 = K.t[jₜ], K.t[jₜ+1]
        tm = t0 + (t1 - t0)/2
        set_t0t1!(ws_D, t0, tm)
        fill_stepMX!(D₁ᵀ, ws_D, rowcomputation, smart_integration)
        set_t0t1!(ws_T, t0, t1)
        fill_stepMX!(Tᵀ, ws_T, rowcomputation, smart_integration)
        set_t0t1!(ws_D, tm, t1)
        fill_stepMX!(D₂ᵀ, ws_D, rowcomputation, smart_integration)
        # Sᵀ = D₁ᵀ Tᵀ D₂ᵀ (serial: the products of sparse matrices are faster serially than in column blocks with threads)
        write_stepMX!(storage_matrix(stepMX[jₜ]), D₁ᵀ * Tᵀ * D₂ᵀ; kwargs...)
    end
end
# The diffusion matrix is block diagonal (sparse), the transport matrix is dense for dense interpolations
function operator_storage(IK)
    n = length(IK.pdf)
    Q_dense = IK.sdestep.steptracer isa TransportStepTracer && get_itp_type(IK.pdf.axes) <: DenseInterpolationType
    Q_dense ? zeros(eltype(IK.pdf.p), n, n) : spzeros(eltype(IK.pdf.p), n, n)
end

# Write P = Sᵀ into the stored matrix A of the step matrix, without the elements |Sᵢⱼ| ≤ max(sparse_tol, sparse_rtol * max|S|)
function write_stepMX!(A::SparseMatrixCSC{Tv,Ti}, P::AbstractMatrix; sparse_tol = 1e-6, sparse_rtol = zero(Tv), kwargs...) where {Tv,Ti}
    P = P isa SparseMatrixCSC ? P : sparse(P)
    vals, rows = nonzeros(P), rowvals(P)
    τ = max(sparse_tol, sparse_rtol * maximum(abs, vals; init = zero(Tv)))
    nnz_A = count(v -> abs(v) > τ, vals)
    nnz_A < typemax(Ti) || throw(ArgumentError("The number of nonzero elements ($nnz_A) is too large for index type $Ti: use SparseMX(index_type = Int64)"))
    colptr, rowval, nzval = SparseArrays.getcolptr(A), rowvals(A), nonzeros(A)
    resize!(rowval, nnz_A)
    resize!(nzval, nnz_A)
    colptr[1] = one(Ti)
    k = 0
    for j in 1:size(P, 2)
        for p in nzrange(P, j)
            if abs(vals[p]) > τ
                k += 1
                rowval[k] = rows[p]
                nzval[k] = vals[p]
            end
        end
        colptr[j+1] = k + 1
    end
    A
end
write_stepMX!(A::Matrix, P::AbstractMatrix; kwargs...) = copyto!(A, P)
function write_stepMX!(A::Matrix, P::SparseMatrixCSC; kwargs...)
    fill!(A, zero(eltype(A)))
    vals, rows = nonzeros(P), rowvals(P)
    for j in 1:size(P, 2), p in nzrange(P, j)
        A[rows[p], j] = vals[p]
    end
    A
end
