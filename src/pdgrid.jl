struct PDGrid{N,k,xeT,xT,pT}
    Δx::xeT
    xs::xT
    p::pT
    # itp::itpT # Chebyshev or equidistant linear grid
end

function PDGrid(sde::SDE{N,k},xs::xsT; Q_equidistant = true) where xsT<:AbstractVector{xT} where xT<:AbstractVector where {N,k}
    if Q_equidistant
        Δx = [x[2]-x[1] for x in xs]
    else
        Δx = nothing
    end
    lens = length.(xs);
    p = zeros(eltype(xs[1]),lens...)
    PDGrid{N,k,typeof(Δx),xsT,typeof(p)}(Δx,xs,p)
end