function compute_stepMX(IK; kwargs...)
    stepMX = initialize_stepMX(eltype(IK.pdf.p), IK.t, length(IK.pdf))

    fill_stepMX_ts!(stepMX, IK; kwargs...)
    stepMX
end

@inline initialize_stepMX(T, ts::AbstractVector{eT}, l) where eT<:Number = [zeros(T,l,l) for _ in 1:(length(ts)-1)]
@inline initialize_stepMX(T, ts::Number, l) = zeros(T, l, l)

function fill_stepMX_ts!(stepMX, IK::IntegrationKernel{kd, sdeT,xT,fT,pdfT, tT}; kwargs...) where {kd, sdeT,xT,fT,pdfT, tT::AbstractArray}
    for jₜ in 1:length(ts)-1
        IK.sdestep.t0[1] = IK.t[jₜ]
        IK.sdestep.t1[1] = IK.t[jₜ+1]
        fill_stepMX!(stepMX, IK)
    end
end
function fill_stepMX_ts!(stepMX, IK::IntegrationKernel{kd, sdeT,xT,fT,pdfT, tT}; kwargs...) where {kd, sdeT,xT,fT,pdfT, tT::Number}
    fill_stepMX!(stepMX, IK)
end

function fill_stepMX!(stepMX, IK::IintegrationKernel{1})
    for (i, idx) in enumerate(IK.temp.idx_it)
        update_state_x1!(IK.sdestep, IK.pdf.axes, idx)
        # TODO
    end
end
function update_state_x1!(sdestep, axes, idx)
    @. sdestep.x1 =  getindex(axes,idx) # TODO: check allocations!
end