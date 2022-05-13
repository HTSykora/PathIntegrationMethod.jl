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
    err = sum(abs, PI.pdf.(_x) .- p_AN(_x)) * ((_x[end] - _x[1])/length(_x));
    err
end

##
euler = Euler()
rk2 = RK2()
rk4 = RK4()

testresults = BitArray(undef,0)
push!(testresults, get_PI_err(51, 0.001, ChebyshevAxis, euler) < 0.0009)
push!(testresults, get_PI_err(51, 0.001, ChebyshevAxis, rk2) < 0.0012)
push!(testresults, get_PI_err(51, 0.001, ChebyshevAxis, rk4) < 0.0012)
reduce(&, testresults)
