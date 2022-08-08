using PathIntegrationMethod

# using QuadGK, Arpack
##
# 1D problem:
f(x,p,t) = x[1]-x[1]^3
g(x,p,t) = sqrt(2)
sde = SDE(f,g)

# Analytic form
_p_AN(x,ε = 1.) = exp(x^2/2 - ε*x^4/4)
function p_AN(xs, ε=1.)
    itg = PathIntegrationMethod.quadgk(x->_p_AN(x,ε) ,xs[1],xs[end])[1]
    _p_AN.(xs,Ref(ε)) ./ itg
end


function get_PI_err(N, Δt, itpf, method; xmin = -3.0, xmax = 3.0, _testN = 10001, _x = LinRange(xmin, xmax, _testN), Tmax = 10.0, kwargs...)
    axisgrid = itpf(_x[1],_x[end],N)
    PI = PathIntegration(sde, method, Δt, axisgrid; kwargs...);
    advance_till_converged!(PI, Tmax = Tmax)
    err1 = sum(abs, PI.pdf.(_x) .- p_AN(_x)) * ((_x[end] - _x[1])/length(_x));
    recompute_PI!(PI)
    advance_till_converged!(PI, Tmax = Tmax)
    err2 = sum(abs, PI.pdf.(_x) .- p_AN(_x)) * ((_x[end] - _x[1])/length(_x));
    err1, err2
end

##
methods = (LinearAxis, CubicAxis, QuinticAxis, ChebyshevAxis, TrigonometricAxis)
err_refs = ([0.0175,0.027,0.027],
            [0.0081,0.011,0.011],
            [0.009,0.0115,0.0115],
            [0.009,0.0115,0.0115],
            [0.009,0.0115,0.0115])

rk1 = Euler()
rk2 = RK2()
rk4 = RK4()

testresults = BitArray(undef,0)
for (i,method) in enumerate(methods)
    err_ref = err_refs[i]
    err_rk1_1, err_rk1_2 = get_PI_err(71, 0.01, method, rk1)
    err_rk2_1, err_rk2_2 = get_PI_err(71, 0.01, method, rk2)
    err_rk4_1, err_rk4_2 = get_PI_err(71, 0.01, method, rk4)

    push!(testresults, isapprox(err_rk1_1, err_rk1_2, atol = 1.5e-8))
    push!(testresults, isapprox(err_rk2_1, err_rk2_2, atol = 1.5e-8))
    push!(testresults, isapprox(err_rk4_1, err_rk4_2, atol = 1.5e-8))
    push!(testresults, err_rk1_1 < err_ref[1])
    push!(testresults, err_rk2_1 < err_ref[2])
    push!(testresults, err_rk4_1 < err_ref[3])
end
reduce(&, testresults)
