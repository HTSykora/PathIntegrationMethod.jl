using Revise, BenchmarkTools
using PathIntegrationMethod

function f2(x,p,t)
    ζ, _ = p # ζ, σ
    -2ζ*x[2] - x[1]
end
function g2(x,p,t)
    p[2] # = σ
end
par = [0.1,0.1]
##
method = Euler()
method = RK2()

W = Wall(x->0.7-0.05x,1.)
W = Wall(0.7,1.)
sde = SDE_VIO(f2,g2,W, par)

@time sdestep = SDEStep(sde, method, 0.1);
step1 = sdestep.sdesteps[1]
step1.x0 .= [0.99, 1.]
step2 = sdestep.sdesteps[2]
step2.x0 .= [0.99, 1.]

@time PathIntegrationMethod.eval_driftstep!(step1)
step1.x1

begin
    step2.xi[1] = W.pos
    step2.xi[2] = step2.x0[2]
    step2.xi2[1] = W.pos
    step2.ti[] = PathIntegrationMethod._Δt(step2)/2

    step2.x1 .= [0.95, -0.5]
    @time PathIntegrationMethod.compute_missing_states_driftstep!(step2, PathIntegrationMethod.update_impact_vio_x!)
end

begin
    @show step2.x0
    @show step2.x1
    @show step2.xi
    @show step2.xi2
end



using Symbolics
W = Wall(x->0.7-0.5x,1.)
W = Wall(0.7,1.)
mdif(W::Wall{<:Number},vi) = zero(vi)
function testSymbolics()
    @variables vi
    @syms r(v)
    expr = Symbolics.derivative(-r(vi)*vi,vi)
    ep = expand_derivatives(substitute(expr,Dict(r(vi)=>W(vi))))
    # return vi,r, expr
    sf = build_function(ep,vi, expression = Val{false})
end
sf = testSymbolics()
vi, r, expr = testSymbolics()

function mdif(W::Wall,_vi)
    @variables vi
    expand_derivatives(Differential(vi)(W(vi)))
end
