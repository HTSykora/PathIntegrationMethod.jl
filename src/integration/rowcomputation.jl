default_rowcomputation() = Threads.nthreads() > 1 ? ThreadedRowComputation() : SerialRowComputation()

n_workspaces(::SerialRowComputation) = 1
n_workspaces(rc::Union{ThreadedRowComputation,BatchRowComputation}) = rc.N_threads

function RowWorkspace(IK, rc::RowComputation, stepMX)
    n = n_workspaces(rc)
    IKs = Vector{typeof(IK)}(undef, n)
    IKs[1] = IK
    for c in 2:n
        IKs[c] = task_copy(IK)
    end
    RowWorkspace(IKs, chunk_ranges(length(IK.pdf.p), n), row_buffers(storage_matrix(stepMX), n))
end
# n blocks of consecutive indices of 1:l with (almost) equal lengths
chunk_ranges(l, n) = [((c-1)*l÷n + 1):(c*l÷n) for c in 1:n]
row_buffers(A::SparseMatrixCSC, n) = (nnz_per_row = zeros(Int, size(A, 2)); [SparseRowBuffer(A, nnz_per_row) for _ in 1:n])
row_buffers(A::AbstractMatrix, n) = nothing

# Where the rows are written: a dense matrix is shared by the chunks, a sparse one is assembled from the buffers
row_outputs(A::AbstractMatrix, ws::RowWorkspace) = fill(A, length(ws.chunks))
function row_outputs(A::SparseMatrixCSC, ws::RowWorkspace)
    for buf in ws.buffers
        empty!(buf.rowval)
        empty!(buf.nzval)
        buf.maxabs[] = zero(buf.maxabs[])
    end
    ws.buffers
end

function set_t0t1!(ws::RowWorkspace, t0, t1)
    for IK in ws.IKs
        set_t0t1!(IK.sdestep, t0, t1)
    end
end

# Compute the rows ws.chunks[c] with ws.IKs[c] for every chunk c
function fill_chunks!(::SerialRowComputation, outs, ws, smart_integration)
    for c in eachindex(ws.chunks)
        fill_rows!(outs[c], ws.IKs[c], ws.chunks[c], smart_integration)
    end
end
# `@threads :static` cannot be nested: inside a threaded region the chunks are computed one after the other
in_threaded_region() = ccall(:jl_in_threaded_region, Cint, ()) != 0
function fill_chunks!(rc::ThreadedRowComputation, outs, ws, smart_integration)
    if Threads.nthreads() == 1 || in_threaded_region()
        return fill_chunks!(SerialRowComputation(), outs, ws, smart_integration)
    end
    with_single_threaded_blas(rc.single_threaded_blas, ws) do
        Threads.@threads :static for c in eachindex(ws.chunks)
            fill_rows!(outs[c], ws.IKs[c], ws.chunks[c], smart_integration)
        end
    end
end
# Polyester passes the arrays of Numbers used in the loop as pointer arrays: everything is passed in a single struct
struct ChunkJob{oT,wT}
    outs::oT
    ws::wT
    smart_integration::Bool
end
fill_chunk!(job::ChunkJob, c) = fill_rows!(job.outs[c], job.ws.IKs[c], job.ws.chunks[c], job.smart_integration)
function fill_chunks!(rc::BatchRowComputation, outs, ws, smart_integration)
    job = ChunkJob(outs, ws, smart_integration)
    with_single_threaded_blas(rc.single_threaded_blas, ws) do
        @batch per=thread for c in eachindex(ws.chunks)
            fill_chunk!(job, c)
        end
    end
end

# Evaluate f() with a single BLAS thread (if Q and the rows use BLAS), then restore the number of BLAS threads
function with_single_threaded_blas(f, Q::Bool, ws)
    n_blas = BLAS.get_num_threads()
    if !Q || n_blas == 1 || !(first(ws.IKs).temp.kernel isa DenseTensorKernel)
        return f()
    end
    BLAS.set_num_threads(1)
    try
        f()
    finally
        BLAS.set_num_threads(n_blas)
    end
end

## Copies of the integration kernel with separate buffers (the data that is only read is shared)
task_copy(x) = deepcopy(x)
function task_copy(IK::IntegrationKernel)
    IntegrationKernel(task_copy(IK.sdestep), IK.f, task_copy(IK.discreteintegrator), IK.t, IK.pdf, task_copy(IK.temp), IK.kwargs)
end

# The SDE (with its parameters) and the compiled Newton steps are shared
function task_copy(s::SDEStep{d,k,m,sdeT,methodT,tracerT,x0T,x1T,tT,tiT,xiT,xi2T}) where {d,k,m,sdeT<:SDE,methodT,tracerT,x0T,x1T,tT,tiT,xiT,xi2T}
    SDEStep{d,k,m,sdeT,methodT,tracerT,x0T,x1T,tT,tiT,xiT,xi2T}(s.sde, deepcopy(s.method), copy(s.x0), copy(s.x1), Ref(s.t0[]), Ref(s.t1[]),
        task_copy(s.steptracer), deepcopy(s.ti), deepcopy(s.xi), deepcopy(s.xi2))
end
task_copy(st::SymbolicNewtonStepTracer) = SymbolicNewtonStepTracer(st.xI_0!, st.x_0!, st.detJI_inv, deepcopy(st.tempI), deepcopy(st.temp))

function task_copy(di::DiscreteIntegrator{dim}) where dim
    DiscreteIntegrator{dim}(copy(di.x), copy(di.w), zero(di.res), zero(di.temp), Ref(di.Q_integrate[]), di.x_ref, di.w_ref, di.rule)
end
task_copy(q::QuadGKIntegrator) = QuadGKIntegrator(copy(q.int_limits), zero(q.res), q.kwargs, Ref(q.Q_integrate[]), zero(q.res0))
function task_copy(di::NonSmoothDiscreteIntegrator{dim,NoDyn,disT}) where {dim,NoDyn,disT}
    NonSmoothDiscreteIntegrator{dim,NoDyn,disT}(map(task_copy, di.discreteintegrators))
end

task_copy(t::IK_temp) = IK_temp(map(zero, t.itpVs), zero(t.itpM), task_copy(t.kernel))
task_copy(::GenericRowKernel) = GenericRowKernel()
task_copy(k::SparseAccumulator) = SparseAccumulator(zero(k.mark), similar(k.touched, 0))
task_copy(k::DenseTensorKernel) = DenseTensorKernel(map(zero, k.Bs), zero(k.c), zero(k.B1c), k.KR isa Nothing ? nothing : zero(k.KR))
