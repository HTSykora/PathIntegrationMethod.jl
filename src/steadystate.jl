"""
    steady_state!(PI; method = :auto, tol = 1e-10, maxiter = 100, near_one_tol = nothing, reset_t = true)

Compute the stationary response PDF of `PI` directly, instead of iterating `advance!` (see `advance_till_converged!`).
The stationary PDF is the eigenvector of the step matrix `S` belonging to the eigenvalue closest to 1.
For time-periodic systems (several step matrices) it is the eigenvector of the period map `S_n⋯S_1`, i.e. the periodic response PDF at the start of the period.

The result is stored in `PI.pdf`, normalised to ``∫p = 1``, and `PI.step_idx` is reset (the next `advance!` uses `S_1`).

Returns `(PI, info)`, where `info = (λ, λ₂, residual, imag_ratio, negative_mass, method)`:
- `λ`: the eigenvalue (≈ 1, the probability kept in a step or period)
- `λ₂`: the eigenvalue second closest to 1 (among the eigenvalues computed by `method`)
- `residual`: ``‖Pp - λp‖/‖p‖``, where `P` is `S` or the period map
- `imag_ratio`: relative size of the imaginary part of the eigenvector (it should be ≈ 0)
- `negative_mass`: ``∫max(-p,0)`` (high-order interpolations can produce slightly negative PDF values)

# Keyword Arguments
- `method = :auto`:
    - `:eigen`: dense eigenvalue decomposition (small problems)
    - `:lu`: sparse LU decomposition of ``S - σI`` (σ ≈ 1) and inverse iteration; time-invariant systems only
    - `:arnoldi`: Krylov–Schur iteration (ArnoldiMethod.jl), only uses matrix–vector products
    - `:auto`: `:eigen` for small problems (at most 64 grid points), `:lu` for a single step matrix, `:arnoldi` otherwise
- `tol = 1e-10`: convergence tolerance of `:lu` and `:arnoldi`
- `maxiter = 100`: maximum number of inverse iterations of `:lu`
- `near_one_tol = nothing`: eigenvalues with ``|λ - 1| ≤`` `near_one_tol` are near 1. If several eigenvalues are near 1,
  the steady state is not well determined (e.g. weakly coupled regions of the state space, such as the wells of a bistable system),
  `advance!` may converge to a different PDF, and a warning is shown.
  The default is ``10|1 - λ|`` (clamped to ``[√eps, 0.1]``): ``|1 - λ|`` is the probability lost or gained in a step,
  i.e. how accurately the discretised step conserves probability, so eigenvalues closer to 1 than a few times this cannot be distinguished reliably from λ.
- `reset_t = true`: set `PI.t` to 0

Warnings are shown if the result may be inaccurate (complex eigenvector, negative PDF values, large residual) or if several eigenvalues are near 1.
"""
function steady_state!(PI::PathIntegration; method::Symbol = :auto, tol = 1e-10, maxiter = 100, near_one_tol = nothing, reset_t = true)
    PI.stepMX isa Nothing && throw(ArgumentError("The step matrix is not computed (pre_compute = false)"))
    Ss = PI.stepMX isa AbstractVector ? PI.stepMX : [PI.stepMX]
    w = vec(quadrature_weights(PI.pdf))

    _method = method === :auto ? auto_steady_state_method(Ss) : method
    # λs: the eigenvalues computed by the method (closest to 1 or largest)
    if _method === :eigen
        λ, v, λs = steady_state_eigen(Ss)
    elseif _method === :lu
        λ, v, λs = steady_state_lu(Ss, vec(PI.pdf.p), w; tol = tol, maxiter = maxiter)
    elseif _method === :arnoldi
        λ, v, λs = steady_state_arnoldi(Ss; tol = tol)
    else
        throw(ArgumentError("Unknown method :$method (use :auto, :eigen, :lu or :arnoldi)"))
    end

    imag_ratio = norm(imag.(v)) / norm(v)
    p = real.(v)
    p ./= dot(w, p) # ∫p = 1, also fixes the sign of the eigenvector
    residual = norm(PeriodMap(Ss) * p - real(λ) * p) / norm(p)
    negative_mass = sum(_w * max(-_p, zero(_p)) for (_w, _p) in zip(w, p))
    if abs(imag(λ)) > 1e-8 || imag_ratio > 1e-6 || negative_mass > 1e-3 || residual > sqrt(tol)
        @warn "The steady-state PDF may be inaccurate" λ imag_ratio negative_mass residual
    end
    λ₂ = second_closest_to_one(λs)
    _near_one_tol = near_one_tol isa Nothing ? clamp(10abs(1 - λ), sqrt(eps(real(float(typeof(λ))))), 0.1) : near_one_tol
    n_near_one = count(μ -> abs(μ - 1) ≤ _near_one_tol, λs)
    if n_near_one ≥ 2
        @warn "Several eigenvalues are near 1: the steady state is not well determined (e.g. weakly coupled regions of the state space), and advance! may converge to a different PDF" λ λ₂ n_near_one near_one_tol = _near_one_tol
    end

    vec(PI.pdf.p) .= p
    PI.step_idx = zero(PI.step_idx)
    if reset_t
        PI.t = zero(PI.t)
    end
    PI, (; λ, λ₂, residual, imag_ratio, negative_mass, method = _method)
end

function auto_steady_state_method(Ss)
    if size(first(Ss), 1) ≤ 64
        :eigen
    elseif length(Ss) == 1
        :lu
    else
        :arnoldi
    end
end

# Eigenvalue closest to 1 and its eigenvector (and all eigenvalues)
function closest_to_one(λs, X)
    i = argmin(abs.(λs .- 1))
    λs[i], X[:, i], λs
end
second_closest_to_one(λs) = length(λs) < 2 ? oftype(float(first(λs)), NaN) : sort(λs, by = μ -> abs(μ - 1))[2]

function steady_state_eigen(Ss)
    P = Matrix(first(Ss))
    for S in Iterators.drop(Ss, 1)
        P = Matrix(S) * P
    end
    E = eigen(P)
    closest_to_one(E.values, E.vectors)
end

# Inverse iteration with the shift σ ≈ 1: converges to the eigenvector of the eigenvalue closest to 1,
# i.e. to the fixed point of `advance!` (S is not exactly probability conserving, so S - I is not singular)
function steady_state_lu(Ss, p0, w; tol = 1e-10, maxiter = 100)
    length(Ss) == 1 || throw(ArgumentError("method = :lu is only available for systems with a single step matrix (use :arnoldi)"))
    S = first(Ss)
    σ = 1 + sqrt(eps(real(eltype(S))))
    F = shifted_lu(S, σ)
    p = iszero(dot(w, p0)) ? ones(eltype(S), length(p0)) : copy(p0)
    p ./= dot(w, p)
    for _ in 1:maxiter
        p_new = F \ p
        p_new ./= dot(w, p_new)
        δ = norm(p_new - p, 1) / norm(p_new, 1)
        p = p_new
        δ < tol && break
    end
    # the two eigenvalues closest to σ (to check whether several eigenvalues are near 1)
    n = size(S, 1)
    decomp, _ = partialschur(ShiftInvert(F, n); nev = min(2, n), which = :LM, tol = sqrt(tol), mindim = min(n, 10), maxdim = min(n, 20))
    λs = σ .+ 1 ./ partialeigen(decomp)[1]
    dot(w, S * p), p, λs # λ = ∫(S p), as ∫p = 1
end
# (S - σI)⁻¹ from the LU decomposition F of S - σI: its largest eigenvalues θ belong to the eigenvalues λ = σ + 1/θ of S closest to σ
struct ShiftInvert{FT}
    F::FT
    n::Int
end
Base.size(A::ShiftInvert) = (A.n, A.n)
Base.size(A::ShiftInvert, i) = A.n
Base.eltype(A::ShiftInvert) = eltype(A.F)
LinearAlgebra.mul!(y::AbstractVector, A::ShiftInvert, x::AbstractVector) = ldiv!(y, A.F, x)
# LU decomposition of S - σI
shifted_lu(S::Transpose{T,<:SparseArrays.AbstractSparseMatrixCSC}, σ) where T = lu(copy(transpose(storage_matrix(S))) - σ*I)
function shifted_lu(S::AbstractMatrix, σ)
    M = Matrix(S) # a new (not transposed) matrix
    M[diagind(M)] .-= σ
    lu!(M)
end

function steady_state_arnoldi(Ss; tol = 1e-10, nev = 4)
    P = PeriodMap(Ss)
    n = size(P, 1)
    # only the eigenpair closest to 1 is needed: its residual is checked in steady_state!
    decomp, _ = partialschur(P; nev = min(nev, n), which = :LM, tol = tol, mindim = min(n, max(10, nev)), maxdim = min(n, max(20, 2nev)))
    closest_to_one(partialeigen(decomp)...)
end

# The period map S_n⋯S_1 of the step matrices Ss as a linear operator
struct PeriodMap{ST,vT}
    Ss::ST
    temp::vT
end
PeriodMap(Ss) = PeriodMap(Ss, zeros(eltype(first(Ss)), size(first(Ss), 1)))
Base.size(P::PeriodMap) = size(first(P.Ss))
Base.size(P::PeriodMap, i) = size(first(P.Ss), i)
Base.eltype(P::PeriodMap) = eltype(first(P.Ss))
function LinearAlgebra.mul!(y::AbstractVector, P::PeriodMap, x::AbstractVector)
    mul!(y, first(P.Ss), x)
    for S in Iterators.drop(P.Ss, 1)
        copyto!(P.temp, y)
        mul!(y, S, P.temp)
    end
    y
end
Base.:*(P::PeriodMap, x::AbstractVector) = mul!(similar(x, eltype(P)), P, x)
