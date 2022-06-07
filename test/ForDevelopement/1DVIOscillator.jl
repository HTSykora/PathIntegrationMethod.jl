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
method = RK2()

W = Wall(x->0.7-0.05x,1.)
W = Wall(0.7,1.)
sde = SDE_VIO(f2,g2,W, par)


@time sdestep = SDEStep(sde, method, 0.1);

sdestep.Q_switch[] = false
step1 = sdestep.sdesteps[1]
step2 = sdestep.sdesteps[2]
step2.x0 .= [0.99, 1.]

begin
    step1.x0 .= [0.99, 1.]
    step1.x0 .= [1., -0.454545454545455]
    step1.x1 .= [0.95, 1.]
    sdestep.Q_switch[] = false
    @time PathIntegrationMethod.compute_missing_states_driftstep!(sdestep)
    @show step1.x0
    @show step1.x1
end

begin
    sdestep.Q_switch[] = true
    sdestep.ID[] = 1

    step2.xi[1] = W.pos
    step2.xi[2] = step2.x0[2]
    step2.xi2[1] = W.pos
    step2.ti[] = PathIntegrationMethod._Δt(step2)/2
    
    step2.x0 .= [0.99, 1.]
    step2.x0 .= [1., 0.64935064935065]

    step2.x1 .= [0.95, -0.5]
    @time PathIntegrationMethod.compute_missing_states_driftstep!(sdestep)

    @show step2.x0
    @show step2.x1
    @show step2.xi
    @show step2.xi2
end;
##
begin
    sdestep.Q_switch[] = true
    sdestep.ID[] = 1

    step2.xi[1] = W.pos
    step2.xi[2] = step2.x0[2]
    step2.xi2[1] = W.pos
    step2.ti[] = PathIntegrationMethod._Δt(step2)/2
    
    step2.x1 .= [0.95, -0.5]
    # step2.x1 .= [0.99, 1.]
    @time PathIntegrationMethod.compute_velocities_to_impact!(step2, sdestep.ID[])

    v_b = step2.steptracer.v_b
    v_a = step2.steptracer.v_a
    @show v_b
    @show v_a
end;

##
# @run sdestep = SDEStep(sde, method, 0.1);
# @run PathIntegrationMethod.compute_velocities_to_impact!(step2, sdestep.ID[])