function fillup!(tpdMX,TT::TransitionTensor{1,1,2})
    zΔt = zero(TT.Δt)
    for (i,x₁) in enumerate(TT.pdgrid.xs[1])
        for (j,x₀) in enumerate(TT.pdgrid.xs[1])
            tpdMX[i,j] = _tp(TT.sde,x₁,TT.Δt, x₀,zΔt)
        end
    end
    tpdMX
end