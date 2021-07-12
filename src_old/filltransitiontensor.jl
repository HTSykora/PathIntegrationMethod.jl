function fillup!(tpdMX, TT::PathIntegrationProblem{1,1,2,T,probT,pdT}) where {T,probT,pdT<:PDGrid{1,1,T2,NTuple{M,xT}}} where {T2,xT<:Axis} where {M}
    temp = zeros(size(tpdMX,2))
    IK = IntegrationKernel(Int[],[0],TT.pdgrid.xs[1],nothing,TT,temp);

    for (idx₁,x₁) in enumerate(TT.pdgrid.xs[1])
        IK.idx₁[1] = idx₁
        quadgk!(IK,IK.temp,IK.xs...)
        tpdMX[idx₁,:] .= IK.temp
        # TODO: intruduce numerical integration scheme
        # tpdMX[idx₁,idx₀] = _tp(TT.sde,x₁,TT.Δt, x₀,zΔt, method = TT.method)
        # tpdMX[idx₁,idx₀] = weight*_tp(TT.sde,x₁,TT.Δt, x₀,zΔt, method = TT.method)
    end

    # for j in 1:size(tpdMX,1)
    #     nrmC = sum(tpdMX[:,j]);
    #     tpdMX[:,j] ./= nrmC
    # end

    tpdMX
end

function fillup!(tpdMX, TT::PathIntegrationProblem{2,2,2})
    zΔt = zero(TT.Δt)
    for (i,x) in enumerate(TT.pdgrid.xs[1])
        for (j,v) in enumerate(TT.pdgrid.xs[2])
            TT.pdgrid.ξ_temp[i,j] = get_ξ(TT.method,TT.sde,TT.Δt,zΔt,x,v)
        end
    end
    # tpdMX .= zero(eltype(tpdMX))
    lx = size(TT.pdgrid,1)
    area = getarea(TT)
    @inbounds for (j₀,v₀) in enumerate(TT.pdgrid.xs[2])
        for (i₀,x₀) in enumerate(TT.pdgrid.xs[1])
            idx₀ = i₀ + (j₀ - 1)*length(TT.pdgrid.xs[1])
            for (j₁,v₁) in enumerate(TT.pdgrid.xs[2])
                for (i₁,x₁) in enumerate(TT.pdgrid.xs[1])
                    idx₁ = i₁ + (j₁ - 1)*lx
                    # tpdMX[idx₁,idx₀] +=_tp(TT.sde,(x₁, v₁),TT.Δt,(x₁-v₀*TT.Δt, v₀), zΔt, method = TT.method)
                    weight = getweight(TT,1,i₀)*getweight(TT,2,j₀)*area
                    tpdMX[idx₁,idx₀] = weight*_tp(TT.sde,(x₁, v₁),TT.Δt,(TT.pdgrid.ξ_temp[i₁,j₀], v₀), zΔt, method = TT.method)
                    # tpdMX[idx₁,idx₀] +=_tp(TT.sde,(x₁, v₁),TT.Δt,(x₀, v₀), zΔt, method = TT.method)
                end
            end
        end
    end
    for j in 1:size(tpdMX,1)
        nrmC = sum(tpdMX[:,j]);
        tpdMX[:,j] ./= nrmC
    end

    tpdMX
end