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

struct WallFunction{mnT, aT, dT, dxT, sfT} <: Function
    mn::mnT
    a::aT
    dx::dxT
    d::dT
    sigmoidfunction::sfT
end
function (w::WallFunction)(x)
    w.mn + w.d * w.sigmoidfunction(-w.a*x + w.dx)
end

function WallFunction(mn, mx, a = 1.; sigmoidfunction::T = tanh) where T<:typeof(tanh)
    d = (mx - mn)/2
    dx = 4.;
    WallFunction(mn + d, a, dx, d, tanh)
end


##
method = Euler()
# method = RK2()

W1 = Wall(0.7,-0.0)
sde1 = SDE_VIO(f2,g2,W1, par)
wf = WallFunction(0.25,0.7,5.0)
W2 = Wall(wf,-0.)
sde2 = SDE_VIO(f2,g2,W2, par)
Δt = 0.01

@time sdestep1 = SDEStep(sde1, method, Δt);
@time sdestep2 = SDEStep(sde2, method, Δt);
i_step1 = sdestep1.sdesteps[2];
i_step2 = sdestep2.sdesteps[2];
##


begin # Check velocity finding
    _x1 = [0.1,0.1]
    i_step1.x1 .= _x1
    i_step2.x1 .= _x1
    @time PathIntegrationMethod.compute_velocities_to_impact!(i_step1)
    @time PathIntegrationMethod.compute_velocities_to_impact!(i_step2)
    v_i1 = i_step1.steptracer.v_i
    v_i2 = i_step2.steptracer.v_i

    println("")
    @show v_i1
    @show v_i2
    println("")
    @show W1(v_i1[1])
    @show W2(v_i2[1])
    
    println("")
    @show PathIntegrationMethod.get_v_beforeimpact(5.,W1)
    @show PathIntegrationMethod.get_v_beforeimpact(5.,W2)
end



begin
    sdestep.ID_dyn[] = 1
    
    step1.x0 .=  [-0.095, -1.]
    # step1.x0 .= [1., -0.454545454545455]
    step1.x1 .= [-0.095, 1.]

    step1.x0 .= [0.,0.0]
    step1.x1 .= [0.,-0.1]
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
    PathIntegrationMethod.set_wall_ID!(step2,1)
    step2.x0 .= [PathIntegrationMethod.get_wall(step2).pos,-0.001]
    step2.x1 .= [PathIntegrationMethod.get_wall(step2).pos, 0.0]
    
    step2.ti[] = Δt
    step2.xi[1] = pos
    step2.xi[2] = step2.x0[2]
    step2.xi2[1] = pos
    step2.xi2[2] = -step2.x0[2]*W(step2.x0[2])
    # step2.ti[] = 0*PathIntegrationMethod._Δt(step2)/2

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

##
# @run sdestep = SDEStep(sde, method, 0.1);
# @run PathIntegrationMethod.compute_velocities_to_impact!(step2, sdestep.ID[])



##
begin
    sdestep.ID_dyn[] = 2
    PathIntegrationMethod.set_wall_ID!(sdestep,1)
    # step2.x0 .= [0.0014430014430014432, -0.014430014430014432]
    # step2.x0 .= rand(2)
    # step2.x1 .= [0.,-0.3]
    # step2.x1 .= [W.pos,-0.000]
    # step2.xi2[2] = step2.x0[2]
    # step2.xi[2] = -step2.xi2[2]/0.7
    # step2.ti[] = 0.00
    
    # step2.xi[1] = W.pos
    # step2.xi[2] = step2.x0[2]
    # step2.xi2[1] = W.pos
    # step2.ti[] = PathIntegrationMethod._Δt(step2)
    
    PathIntegrationMethod.update_dyn_state_xs!(step2,[1e-4,0.])
    @show step2.x0
    @show step2.x1
    @show step2.xi
    @show step2.xi2
    @show step2.ti[]

    @time PathIntegrationMethod.compute_initial_states_driftstep!(step2)
    @show step2.x0
    @show step2.x1
    @show step2.xi
    @show step2.xi2
    @show step2.ti[]
    @show step2.steptracer.temp

    PathIntegrationMethod.update_dyn_state_xs!(step2,[1e-4,0.0])
    @time PathIntegrationMethod.compute_missing_states_driftstep!(step2)
    @show step2.x0
    @show step2.x1
    @show step2.xi
    @show step2.xi2
    @show step2.ti[]
    @show step2.steptracer.tempI
end;
