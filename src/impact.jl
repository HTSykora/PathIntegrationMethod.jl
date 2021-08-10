function walls(_d1,_d2)
    d1, d2 = minmax(_d1, _d2)
    Wall(d1,Int8(-1)), Wall(d2,Int8(1))
end
function walls(d)
    @assert d>0 "d<0!"
    walls(-d,d)
end

function hit(x,w::wT) wT<:Tuple{wT1,wT2} where {wT1<:Wall, wT2<:Wall}
    if x<w[1].d
        return true, 1
    elseif x>w[2].d
        return true, 2
    else
        return false, 0
    end
end