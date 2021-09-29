function Quadrature(N) 
    # x,w = gausslegendre(N)
    x,w = gausslobatto(N)
    Quadrature(x,w)
end

function (q::Quadrature)(f, vals)
    for (i,x) in enumerate(q.x)
        vals[i] = f(x)
    end
    sum(fx*w for (fx,w) in zip(vals,q.w))
end

