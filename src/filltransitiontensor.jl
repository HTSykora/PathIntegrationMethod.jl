function fillup!(tpdMX, TT::TransitionTensor{1,1,2})
    zΔt = zero(TT.Δt)
    for (j,x₀) in enumerate(TT.pdgrid.xs[1])
        for (i,x₁) in enumerate(TT.pdgrid.xs[1])
            tpdMX[i,j] = _tp(TT.sde,x₁,TT.Δt, x₀,zΔt, method = TT.method)
        end
    end
    tpdMX
end

function fillup!(tpdMX, TT::TransitionTensor{2,2,2})
    zΔt = zero(TT.Δt)
    for (j,v) in enumerate(TT.pdgrid.xs[2])
        for (i,x) in enumerate(TT.pdgrid.xs[1])
            TT.pdgrid.ξ_temp[i,j] = get_ξ(TT.method,TT.sde,TT.Δt,zΔt,x,v)
        end
    end
    lx = size(TT.pdgrid,1)
    for (j₁,v₁) in enumerate(TT.pdgrid.xs[2])
        for (i₁,x₁) in enumerate(TT.pdgrid.xs[1])
            for (j₀,v₀) in enumerate(TT.pdgrid.xs[2])
                for (i₀,x₀) in enumerate(TT.pdgrid.xs[1])
                    i = i₀ + (j₀ - 1)*lx
                    j = i₁ + (j₁ - 1)*lx
                    tpdMX[i,j] = _tp(TT.sde,(TT.pdgrid.ξ_temp[i₀,j₀], v₀),TT.Δt,(x₁, v₁), zΔt, method = TT.method)
                end
            end
        end
    end
    tpdMX
end