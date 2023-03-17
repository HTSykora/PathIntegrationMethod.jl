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
method = Euler()
# method = RK2()

W = Wall(x->0.7-0.05x,-1.)
W = Wall(0.7,-0.0)
sde = SDE_VIO(f2,g2,W, par)
Δt = 0.1
axisgrid = (QuinticAxis(W.pos,1.,101), QuinticAxis(-3.0,3.,51))

discreteintegrator = (GaussLegendreIntegrator(31),GaussLegendreIntegrator(201))
@time IK = PathIntegration(sde, method, Δt, axisgrid..., discreteintegrator = discreteintegrator, extract_IK = Val{true}());
step1, step2 = IK.sdestep.sdesteps
di1, di2 = IK.discreteintegrator.discreteintegrators
##

idx = (1,28);
getindex.(axisgrid,idx) |>println

PathIntegrationMethod.update_IK_state_x1_by_idx!(IK, idx)
PathIntegrationMethod.update_dyn_state_xs!(IK)


PathIntegrationMethod.rescale_discreteintegrator!(IK, int_limit_thickness_multiplier = 10, impact_int_limit_thickness_multiplier = 10)

begin
    @show IK.x1
    for (i,di) in enumerate(IK.discreteintegrator.discreteintegrators)
        if di.Q_integrate[]
            println("$i = ($(round(di.x[1],digits=5)), $(round(di.x[end],digits=5)))") 
        end
    end
end
begin
    println("1 = ($(round(di1.x[1],digits=5)), $(round(di1.x[end],digits=5)))") 
    @show IK.x1
    step1.x1[1] = IK.x1[1]
    step1.x0[2] = di1.x[1]
    @time PathIntegrationMethod.compute_missing_states_driftstep!(step1)
    @show step1.x0
    @show step1.x1
    @show step1.steptracer.tempI
    @show PathIntegrationMethod.get_detJinv(IK.sdestep[1])
end
begin
    @show IK.x1
    step2.x1[1] = IK.x1[1]
    step2.x0[2] = di2.x[end]
    @time PathIntegrationMethod.compute_missing_states_driftstep!(step2)
    @show step2.x0
    @show step2.x1
    @show step2.xi
    @show step2.xi2
    @show step2.ti[]
    @show step2.steptracer.tempI
    @show PathIntegrationMethod.get_detJinv(IK.sdestep[2])
end


PathIntegrationMethod.rescale_discreteintegrator!(IK, int_limit_thickness_multiplier = 6, impact_int_limit_thickness_multiplier = 999)
PathIntegrationMethod.get_IK_weights!(IK; IK.kwargs...)
sum(abs,IK.temp.itpM) |> println

