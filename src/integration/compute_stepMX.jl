function SparseMX(;threaded = true, sparse_tol::T = 1e-6, kwargs...) where T
    SparseMX{threaded,T}(threaded,sparse_tol)
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
@inline initialize_stepMX(T::DataType, l::Integer, ::SparseMX) = spzeros(T, l, l)
@inline initialize_stepMX(T::DataType, l::Integer, ::DenseMX) = zeros(T, l, l)

function fill_stepMX_ts!(stepMX::AbstractVector{aT}, IK::IntegrationKernel{kd, sdeT,x1T, diT,fT,pdfT, tT}; kwargs...) where {kd, sdeT,x1T, diT,fT,pdfT, aT<:AbstractMatrix{T},tT<:AbstractArray} where T<:Number
    for jₜ in 1:length(IK.t)-1
        set_t0t1!(IK.sdestep, IK.t[jₜ], IK.t[jₜ+1])
        fill_stepMX!(stepMX[jₜ], IK; kwargs...)
    end
end
function fill_stepMX_ts!(stepMX, IK::IntegrationKernel{kd, sdeT,x1T, diT,fT,pdfT, tT}; kwargs...) where {kd, sdeT,x1T, diT,fT,pdfT, tT<:Number}
    fill_stepMX!(stepMX, IK; kwargs...)
end

fill_stepMX!(stepMX::Transpose, IK; kwargs...) = fill_stepMX!(stepMX.parent, IK; kwargs...)
function fill_stepMX!(stepMX, IK; smart_integration = true, kwargs...)
    for (i, idx) in enumerate(dense_idx_it(IK))
        update_IK_state_x1!(IK, idx)
        update_dyn_state_x1!(IK, idx)
        if smart_integration
            rescale_discreteintegrator!(IK; IK.kwargs...)
        end
        get_IK_weights!(IK; IK.kwargs...)
        fill_to_stepMX!(stepMX,IK,i; IK.kwargs...)
    end
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
@inline function fill_to_stepMX!(stepMX::AbstractSparseMatrix,IK,i; sparse_tol = 1e-6, kwargs...)
    for (j,val) in enumerate(IK.temp.itpM)
        if abs(val) > sparse_tol
            # stepMX[i,j] = val
            stepMX[j,i] = val
        end
        # ? fill by rows and multiply from the right when advancing time
    end
    nothing
end

make_sparse(stepMX::AbstractVector{T}) where T<:Number = sparse(stepMX)
make_sparse(stepMX::AbstractVector{T}) where T<:AbstractArray = sparse.(stepMX)

function rescale_discreteintegrator!(IK::IntegrationKernel{1,dyn}; kwargs...) where dyn <:SDEStep{d,k,m} where {d,k,m}
    compute_initial_states_driftstep!(IK.sdestep; IK.kwargs...)
    rescale_discreteintegrator!(IK.discreteintegrator, IK.sdestep, IK.pdf; kwargs...)
end