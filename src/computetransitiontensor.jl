function computeintegrationmatrix(sde::AbstractSDE{N,k},pdgrid,Δt,method; kwargs...) where {N,k}
    warn("Pre-computation of the transitional matrices of an $(N)-dimensional problem with $(N-k+1) noises is not available yet, thus `nothing` is returned!")
    return nothing
end

function computeintegrationmatrix(sde::AbstractSDE{1,1},pdgrid::PDGrid{1,1,T},Δt,method; kwargs...) where {T}
    tpdMX = zeros(T,length(pdgrid.xs...),length(pdgrid.xs...))
    temp = similar(tpdMX, size(tpdMX,2))
    IK = IntegrationKernel(sde,nothing,pdgrid.xs[1],nothing,[0],pdgrid,zero(Δt),Δt,method,temp);
    for idx₁ in eachindex(pdgrid.xs[1])
        IK.idx₁[1] = idx₁
        get_IK_weights!(IK)
        tpdMX[idx₁,:] .= IK.temp
    end

    # for j in 1:size(tpdMX,1)
    #     nrmC = sum(tpdMX[:,j]);
    #     tpdMX[:,j] ./= nrmC
    # end

    return tpdMX
end
