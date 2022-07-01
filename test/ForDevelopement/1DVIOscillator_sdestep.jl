using Revise, BenchmarkTools
using PathIntegrationMethod

using BenchmarkTools

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
# method = RK2()

W = Wall(x->0.7-0.05x,-1.)
W = Wall(0.7,-0.1)
sde = SDE_VIO(f2,g2,W, par)


@time sdestep = SDEStep(sde, method, 0.1);
step1 = sdestep.sdesteps[1]
step2 = sdestep.sdesteps[2]

##

begin
    sdestep.ID_dyn[] = 1

    step1.x0 .=  [-0.095, -1.]
    # step1.x0 .= [1., -0.454545454545455]
    step1.x1 .= [-0.095, 1.]
    @time PathIntegrationMethod.compute_missing_states_driftstep!(step1)
    @show step1.x0
    @show step1.x1
end

begin
    step1.x1[2] += 0.1
    PathIntegrationMethod._getTPDF(sdestep,step1.x0[2],step1.x1)
end
@run PathIntegrationMethod._getTPDF(sdestep,step1.x0[2],step1.x1)

begin
    sdestep.ID_dyn[] = 2
    sdestep.ID_aux[] = 1
    step2.x0 .= [-0.05, -1.]
    step2.x1 .= [-0.075, -0.714285714285715]



    # step2.xi[1] = W.pos
    # step2.xi[2] = step2.x0[2]
    # step2.xi2[1] = W.pos
    # step2.ti[] = PathIntegrationMethod._Δt(step2)
    
    @time PathIntegrationMethod.compute_missing_states_driftstep!(step2)

    @show step2.x0
    @show step2.x1
    @show step2.xi
    @show step2.xi2
    @show step2.ti[]
end;

begin
    PathIntegrationMethod._getTPDF(sdestep,step2.x0[2],step2.x1 +[0,0.1])
end
@run PathIntegrationMethod._getTPDF(sdestep,step2.x0[2],step2.x1 +[0,0.1])
##
begin
    sdestep.ID_aux[] = 1
    step2.x0 .= [-0.98, -1.]
    step2.x1 .= [-0.099, -3.]

    step2.xi[1] = W.pos
    step2.xi[2] = step2.x0[2]
    step2.xi2[1] = W.pos
    step2.ti[] = PathIntegrationMethod._Δt(step2)/2
    
    # step2.x1 .= [0.95, 0.5]
    # step2.x1 .= [0.99, 1.]
    @time PathIntegrationMethod.compute_velocities_to_impact!(step2, sdestep.ID_aux[])

    v_i = step2.steptracer.v_i
    @show v_i # [v_beforeimpact, v_afterimpact, v1]
end;

##
# @run sdestep = SDEStep(sde, method, 0.1);
# @run PathIntegrationMethod.compute_velocities_to_impact!(step2, sdestep.ID[])


PI = PathIntegration(sde, method, Δt, axisgrid; kwargs...);