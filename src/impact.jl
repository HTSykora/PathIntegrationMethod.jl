function Wall(_r,d,id)
    r = Scalar_Or_Function(_r)
    Wall{typeof(r),typeof(d),typeof(id)}(r,d,id)
end

function walls(_d1,_d2, _r1, _r2)
    d1, d2 = minmax(_d1, _d2)
    if d1 == _d1
        r1, r2 = _r1, _r2
    else
        r1, r2 = _r2, _r1
    end
    Wall(r1,d1,Int8(-1)), Wall(r2,d2,Int8(1))
end
function symmetricwalls(d,r)
    @assert d>0 "d<0!"
    walls(-d, d, r, r)
end

function Q_hit(x,w::wT) wT<:Tuple{wT1,wT2} where {wT1<:Wall, wT2<:Wall}
    if x<w[1].d
        return true, 1
    elseif x>w[2].d
        return true, 2
    else
        return false, 0
    end
end

## VI problem initialize
function create_symmetric_VI_PDGrid(sde::SDE_Oscillator1D, v_ax::aT, d, r, Nₓ::Integer; x_interpolation = :chebyshev, kwargs...) where aT<:Axis
    walls = symmetricwalls(d,r);
    x_ax = Axis(-d, d, Nₓ, interpolation = x_interpolation)
    vi_sde = SDE_VI_Oscillator1D(sde,walls)
    PDGrid(vi_sde, x_ax, v_ax; kwargs...)
end