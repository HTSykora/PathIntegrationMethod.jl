using Revise
using PathIntegrationMethod
using PyPlot; pygui(true)
# using QuadGK, Arpack
##
# 1D problem:
f(x,p,t) = x[1]-x[1]^3
g(x,p,t) = 1. #sqrt(2)
# f(x,p,t) = 0. #x[1]-x[1]^3
# g(x,p,t) = 1.
sde = SDE(f,g)


method = RK2()
Δt = 0.001;
sdestep = SDEStep(sde, method, Δt) 
axis = QuinticAxis(-3.5,2.2,101)
# condition = DoubleBarrierCondition(1,2.)
condition = SingleBarrierCondition(1,1.0)
T_mx = MeanFirstPassageTime(sdestep, condition, axis, di_N = 31)

# begin
#     Tf = InterpolatedFunction(axis)
#     Tf.p[T_mx[1]] .= (I-T_mx[2])\fill(Δt,size(T_mx[2],1))
# end

begin
    xs = LinRange_fromaxis(axis,201)
    # figure(1); clf()
    plot(xs, T_mx.T.(xs))
    ylim(bottom=0)#,top=4.05)
end


# Analytic form of first hitting time for Wiener process
# _p_AN(x,ε = 1.) = exp(x^2/2 - ε*x^4/4)
# function p_AN(xs, ε=1.)
#     itg = PathIntegrationMethod.quadgk(x->_p_AN(x,ε) ,xs[1],xs[end])[1]
#     _p_AN.(xs,Ref(ε)) ./ itg
# end
# begin
#     xs = LinRange_fromaxis(axis,201)
#     figure(2); clf()
#     plot(xs,p_AN(xs, 1.))
# end
##


# begin
#     _myf= x->exp(-x^2)
#     myf = InterpolatedFunction(axis,f=_myf)
#     xs = LinRange(axis[1]-1,axis[end]+1,101)
#     figure(1); clf()
#     plot(xs, myf.(xs, allow_extrapolation = true))
#     plot(xs, _myf.(xs))
# end