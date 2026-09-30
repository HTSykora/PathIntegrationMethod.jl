function SparseMX(;threaded = true, sparse_tol = 1e-6, sparse_rtol = 0.0, index_type::Type{Ti} = Int, kwargs...) where Ti<:Integer
    tol, rtol = promote(float(sparse_tol), float(sparse_rtol))
    SparseMX{threaded,typeof(tol),Ti}(threaded,tol,rtol)
end
function get_stepMXtype(sde::AbstractSDE{d},::T; multithreaded_sparse = true, kwargs...) where {d,T}
    SparseMX(; threaded = multithreaded_sparse, kwargs...)
end
function get_stepMXtype(sde::sdeT,::Val{SparseInterpolationType}; multithreaded_sparse = true, kwargs...) where sdeT <: Union{AbstractSDE{1},AbstractSDE{2}}
    SparseMX(; threaded = multithreaded_sparse, kwargs...)
end
function get_stepMXtype(sde::sdeT,::Val{DenseInterpolationType}; kwargs...) where sdeT <: Union{AbstractSDE{1},AbstractSDE{2}}
    DenseMX()
end
get_tol(::DenseMX) = zero(Float64)
get_tol(m::SparseMX) = m.tol
get_rtol(::DenseMX) = zero(Float64)
get_rtol(m::SparseMX) = m.rtol


function compute_stepMX(IK; stepMXtype = DenseMX(), kwargs...)
    stepMX = initialize_stepMX(eltype(IK.pdf.p), IK.t, length(IK.pdf),stepMXtype)

    fill_stepMX_ts!(stepMX, IK; kwargs...)
    get_final_stepMX_form(stepMX, stepMXtype)
    # stepMX
end

@inline get_final_stepMX_form(stepMX::Union{AbstractMatrix{T},AbstractVector{aT}}, ::DenseMX) where aT<:AbstractMatrix{T} where T<:Number = stepMX
@inline get_final_stepMX_form(stepMX::AbstractVector{aT}, mts::SparseMX) where aT<:AbstractSparseMatrix{T} where T<:Number = get_final_stepMX_form.(stepMX,Ref(mts))
@inline function get_final_stepMX_form(stepMX::AbstractSparseMatrix{T},::SparseMX{true}) where T<:Number
    transpose(ThreadedSparseMatrixCSC(stepMX))
end
@inline function get_final_stepMX_form(stepMX::AbstractSparseMatrix{T},::SparseMX{false}) where T<:Number
    transpose(stepMX)
end

@inline initialize_stepMX(T, ts::AbstractVector{eT}, l::Integer, stepMXtype) where eT<:Number = [initialize_stepMX(T,l,stepMXtype) for _ in 1:(length(ts)-1)]
# @inline initialize_stepMX(T, ts::Number, l::Integer, stepMXtype) = [initialize_stepMX(T, l, stepMXtype)]
@inline function initialize_stepMX(T::DataType, l::Integer, ::SparseMX{tf,tolT,Ti}) where {tf,tolT,Ti}
    l < typemax(Ti) || throw(ArgumentError("The number of grid points ($l) is too large for index type $Ti"))
    spzeros(T, Ti, l, l)
end
@inline initialize_stepMX(T::DataType, l::Integer, ::DenseMX) = zeros(T, l, l)

function fill_stepMX_ts!(stepMX::AbstractVector{aT}, IK::IntegrationKernel{kd, sdeT,x1T, diT,fT,pdfT, tT}; rowcomputation = SerialRowComputation(), smart_integration = true, kwargs...) where {kd, sdeT,x1T, diT,fT,pdfT, aT<:AbstractMatrix{T},tT<:AbstractArray} where T<:Number
    ws = RowWorkspace(IK, rowcomputation, first(stepMX))
    for jₜ in 1:length(IK.t)-1
        set_t0t1!(ws, IK.t[jₜ], IK.t[jₜ+1])
        fill_stepMX!(stepMX[jₜ], ws, rowcomputation, smart_integration)
    end
end
function fill_stepMX_ts!(stepMX, IK::IntegrationKernel{kd, sdeT,x1T, diT,fT,pdfT, tT}; rowcomputation = SerialRowComputation(), smart_integration = true, kwargs...) where {kd, sdeT,x1T, diT,fT,pdfT, tT<:Number}
    fill_stepMX!(stepMX, RowWorkspace(IK, rowcomputation, stepMX), rowcomputation, smart_integration)
end

# The matrix where the rows are stored: the rows of S are the columns of the stored (CSC) matrix Sᵀ of a sparse step matrix
storage_matrix(stepMX::Transpose) = storage_matrix(parent(stepMX))
storage_matrix(stepMX::ThreadedSparseMatrixCSC) = stepMX.A
storage_matrix(stepMX::AbstractMatrix) = stepMX

function fill_stepMX!(stepMX, ws::RowWorkspace, rc::RowComputation, smart_integration)
    A = storage_matrix(stepMX)
    fill_chunks!(rc, row_outputs(A, ws), ws, smart_integration)
    finalise_stepMX!(A, ws)
end
finalise_stepMX!(A::AbstractMatrix, ws) = A
finalise_stepMX!(A::SparseMatrixCSC, ws) = assemble_stepMX!(A, ws.buffers; first(ws.IKs).kwargs...)

function fill_rows!(out, IK, rows, smart_integration)
    CI = CartesianIndices(IK.pdf.p)
    for i in rows
        compute_row!(out, IK, i, CI[i], smart_integration)
    end
    nothing
end
function compute_row!(out, IK, i, idx, smart_integration)
    update_IK_state_x1!(IK, idx)
    update_dyn_state_x1!(IK, idx)
    if smart_integration
        rescale_discreteintegrator!(IK; IK.kwargs...)
    end
    get_IK_weights!(IK; IK.kwargs...)
    fill_to_stepMX!(out,IK,i; IK.kwargs...)
end

function update_IK_state_x1!(IK::IntegrationKernel{kd,dyn}, idx) where dyn <:SDEStep{d,k,m} where {kd,d,k,m}
    for i in 1:d
        IK.x1[i] = getindex(IK.pdf.axes[i],idx[i])
    end
end

function update_dyn_state_x1!(IK::IntegrationKernel{kd,dyn}, idx) where dyn <:SDEStep{d,k,m} where {kd, d,k,m}

    update_dyn_state_x1!(IK.sdestep,IK.x1)
    # IK.sdestep.x1 .=  getindex.(IK.pdf.axes,idx) # ? check allocations
    # for i in 1:d
    #     IK.sdestep.x1[i] = getindex(IK.pdf.axes[i],idx[i])
    # end
end
function update_dyn_state_x1!(sdestep::SDEStep{d,k,m}, x1) where {d,k,m}
    sdestep.x1 .= x1
    sdestep.x0 .= x1
end

@inline function fill_to_stepMX!(stepMX::AbstractMatrix,IK,i; kwargs...)
    for j in eachindex(IK.temp.itpM)
        stepMX[i,j] = IK.temp.itpM[j]
        # ? fill by rows and multiply from the right when advancing time
    end
    nothing
end

# Rows of the step matrix in CSC order: row i of S (column i of Sᵀ) has the elements
# rowval[k], nzval[k] for k in (sum(nnz_per_row[1:i-1]) + 1):sum(nnz_per_row[1:i])
struct SparseRowBuffer{Ti,Tv}
    rowval::Vector{Ti}
    nzval::Vector{Tv}
    nnz_per_row::Vector{Int}
    maxabs::Base.RefValue{Tv} # max|S_ij| of the stored elements
end
SparseRowBuffer(A::SparseMatrixCSC{Tv,Ti}, nnz_per_row = zeros(Int, size(A, 2))) where {Tv,Ti} = SparseRowBuffer(Ti[], Tv[], nnz_per_row, Ref(zero(Tv)))

@inline function fill_to_stepMX!(buf::SparseRowBuffer,IK,i; sparse_tol = 1e-6, kwargs...)
    n0 = length(buf.rowval)
    push_row!(buf, IK.temp.kernel, IK.temp.itpM, sparse_tol)
    buf.nnz_per_row[i] = length(buf.rowval) - n0
    nothing
end
function push_row!(buf::SparseRowBuffer, ::AbstractRowKernel, itpM, tol)
    for (j,val) in enumerate(itpM)
        push_element!(buf, j, val, tol)
    end
end
function push_row!(buf::SparseRowBuffer, spa::SparseAccumulator, itpM, tol)
    sort!(spa.touched)
    for j in spa.touched
        push_element!(buf, j, itpM[j], tol)
    end
end
@inline function push_element!(buf::SparseRowBuffer, j, val, tol)
    if abs(val) > tol
        push!(buf.rowval, j)
        push!(buf.nzval, val)
        buf.maxabs[] = max(buf.maxabs[], abs(val))
    end
    nothing
end

# Write the rows stored in `bufs` (consecutive blocks of rows, in order) into the CSC matrix `A` (= Sᵀ)
# Elements with |S_ij| ≤ max(sparse_tol, sparse_rtol * max|S|) are dropped
function assemble_stepMX!(A::SparseMatrixCSC{Tv,Ti}, bufs; sparse_tol = 1e-6, sparse_rtol = zero(Tv), kwargs...) where {Tv,Ti}
    τ = max(sparse_tol, sparse_rtol * maximum(buf -> buf.maxabs[], bufs))
    Q_filter = τ > sparse_tol # the buffers only contain elements with |S_ij| > sparse_tol
    nnz_max = sum(buf -> length(buf.rowval), bufs)
    nnz_max < typemax(Ti) || throw(ArgumentError("The number of nonzero elements ($nnz_max) is too large for index type $Ti: use SparseMX(index_type = Int64)"))
    
    colptr, rowval, nzval = SparseArrays.getcolptr(A), rowvals(A), nonzeros(A)
    resize!(rowval, nnz_max)
    resize!(nzval, nnz_max)
    nnz_per_row = first(bufs).nnz_per_row
    colptr[1] = one(Ti)
    k = 0 # number of elements written
    b, pos = 1, 0 # buffer and position of the last element read
    for i in eachindex(nnz_per_row)
        m = nnz_per_row[i]
        if m > 0
            while pos + m > length(bufs[b].rowval) # row i is in a later buffer
                b, pos = b + 1, 0
            end
            buf = bufs[b]
            for p in pos+1:pos+m
                if !Q_filter || abs(buf.nzval[p]) > τ
                    k += 1
                    rowval[k] = buf.rowval[p]
                    nzval[k] = buf.nzval[p]
                end
            end
            pos += m
        end
        colptr[i+1] = k + 1
    end
    resize!(rowval, k)
    resize!(nzval, k)
    A
end

make_sparse(stepMX::AbstractVector{T}) where T<:Number = sparse(stepMX)
make_sparse(stepMX::AbstractVector{T}) where T<:AbstractArray = sparse.(stepMX)

function rescale_discreteintegrator!(IK::IntegrationKernel{1,dyn}; kwargs...) where dyn <:SDEStep{d,k,m} where {d,k,m}
    compute_initial_states_driftstep!(IK.sdestep; IK.kwargs...)
    rescale_discreteintegrator!(IK.discreteintegrator, IK.sdestep, IK.pdf; kwargs...)
end