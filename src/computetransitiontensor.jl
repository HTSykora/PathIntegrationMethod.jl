function computeintegrationmatrix(sde::AbstractSDE{N,k},pdgrid,Δt,method; kwargs...) where {N,k}
    warn("Pre-computation of the transitional matrices of an $(N)-dimensional problem with $(N-k+1) noises is not available yet, thus `nothing` is returned!")
    return nothing
end

function computeintegrationmatrix(sde::AbstractSDE{1,1},pdgrid::PDGrid{1,1,T},Δt,method; kwargs...) where {T}
    tpdMX = zeros(T,length(pdgrid),length(pdgrid))
    temp = similar(pdgrid)
    IK = IntegrationKernel(sde,nothing,pdgrid.xs[1],nothing,[0],pdgrid,zero(Δt),Δt,method,temp);
    for idx₁ in eachindex(pdgrid.xs[1])
        IK.idx₁[1] = idx₁
        get_IK_weights!(IK)
        tpdMX[idx₁,:] .= IK.temp
    end

    return tpdMX
end

function computeintegrationmatrix(sde::SDE_Oscillator1D,pdgrid::PDGrid{2,2,T},Δt,method; kwargs...) where {T}
    # Prepare IK
    tpdMX = zeros(T,length(pdgrid),length(pdgrid))
    temp = similar(pdgrid)
    IK = IntegrationKernel(sde,nothing,pdgrid.xs[2],nothing,[0,0],pdgrid,zero(Δt),Δt,method,temp);
    # Fill the matrix representation of the transition tensor (tpdMX)
    idx_it = Base.Iterators.product(eachindex.(pdgrid.xs)...)

    for (i,idx₁) in enumerate(idx_it)
        IK.idx₁ .= idx₁
        get_IK_weights!(IK)
        for j in eachindex(IK.temp)
            tpdMX[i,j] = IK.temp[j] # rework it to use the transpose!
        end
    end

    # for i in 1:length(pdgrid)
    #     _norm = sum(view(tpdMX,:,i))
    #     for j in 1:length(pdgrid)
    #         tpdMX[j,i] = tpdMX[j,i] / _norm
    #     end
    # end
    return tpdMX
end
# function computeintegrationmatrix(sde::SDE_Oscillator1D,pdgrid::PDGrid{2,2,T},Δt,method; kwargs...) where {T}
#     # Prepare IK
#     tpdMX = zeros(T,size(pdgrid)...,size(pdgrid)...)
#     temp = similar(pdgrid)
#     IK = IntegrationKernel(sde,nothing,pdgrid.xs[2],nothing,[0,0],pdgrid,zero(Δt),Δt,method,temp);

#     # Fill the matrix representation of the transition tensor (tpdMX)
#     idx_it = Base.Iterators.product(eachindex.(pdgrid.xs)...)

#     for (i,idx₁) in enumerate(idx_it)
#         IK.idx₁ .= idx₁
#         get_IK_weights!(IK)
#         tpdMX[:,:,idx₁...] .= IK.temp # rework it to use the transpose!
#     end

#     return tpdMX
# end
