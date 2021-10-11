function compute_stepMX(IK; kwargs...)
    stepMX = initialize_stepMX(eltype(IK.pdf.p), IK.t, length(IK.pdf))

    fill_stepMX_ts!(stepMX, IK; kwargs...)
    stepMX
end

@inline initialize_stepMX(T, ts::AbstractVector{eT}, l) where eT<:Number = [zeros(T,l,l) for _ in 1:(length(ts)-1)]
@inline initialize_stepMX(T, ts::Number, l) = zeros(T, l, l)

function fill_stepMX_ts!(stepMX, IK::IntegrationKernel{kd, sdeT,x1T, diT,fT,pdfT, tT}; kwargs...) where {kd, sdeT,x1T, diT,fT,pdfT, tT<:AbstractArray}
    for jₜ in 1:length(IK.ts)-1
        IK.sdestep.t0[1] = IK.t[jₜ]
        IK.sdestep.t1[1] = IK.t[jₜ+1]
        fill_stepMX!(stepMX, IK)
    end
end
function fill_stepMX_ts!(stepMX, IK::IntegrationKernel{kd, sdeT,x1T, diT,fT,pdfT, tT}; kwargs...) where {kd, sdeT,x1T, diT,fT,pdfT, tT<:Number}
    fill_stepMX!(stepMX, IK)
end

function fill_stepMX!(stepMX, IK::IntegrationKernel)
    for (i, idx) in enumerate(dense_idx_it(IK))
        update_IK_state_x1!(IK, idx)
        update_dyn_state_x1!(IK, idx)
        rescale_discreteintegrator!(IK)
        get_IK_weights!(IK)
        fill_to_stepMX!(stepMX,IK,i)
    end
end

function update_IK_state_x1!(IK::IntegrationKernel{kd,dyn}, idx) where dyn <:SDEStep{d,k,m} where {kd,d,k,m}
    for i in 1:d
        IK.x1[i] = getindex(IK.pdf.axes[i],idx[i])
    end
end

function update_dyn_state_x1!(IK::IntegrationKernel{kd,dyn}, idx) where dyn <:SDEStep{d,k,m} where {kd, d,k,m}

    IK.sdestep.x1 .= IK.x1
    # IK.sdestep.x1 .=  getindex.(IK.pdf.axes,idx) # ? check allocations
    # for i in 1:d
    #     IK.sdestep.x1[i] = getindex(IK.pdf.axes[i],idx[i])
    # end
end

@inline function fill_to_stepMX!(stepMX,IK,i)
    for j in eachindex(IK.temp.itpM)
        stepMX[i,j] = IK.temp.itpM[j]
        # ? fill by rows and multiply from the right when advancing time
    end
    nothing
end

function rescale_discreteintegrator!(IK::IntegrationKernel{1,dyn}; σ_mult = 10, kwargs...) where dyn <:SDEStep{d,k,m} where {kd,d,k,m}
    compute_initial_states_driftstep!(IK.sdestep)
    σ = sqrt(_Δt(IK.sdestep))*abs(IK.sdestep.sde.g(d, IK.sdestep.x0,_par(IK.sdestep),_t0(IK.sdestep)))
    mn = max(IK.pdf.axes[d][1],IK.sdestep.x0[d] - σ_mult*σ)
    mx = min(IK.pdf.axes[d][end],IK.sdestep.x0[d] + σ_mult*σ)

    rescale_to_limits!(IK.discreteintegrator,mn,mx)
end