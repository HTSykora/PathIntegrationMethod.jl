using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];
##

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
par = [0.1,0.5]
##
# method = RK2()

#  W = Wall(x->0.7-0.001x^2,-0.)
W = Wall(0.7,-0.0)
sde = SDE_VIO(f2,g2,W, par)
axisgrid = (QuinticAxis(W.pos,4.,101), QuinticAxis(-3.0,3.,101))

# Ws = (Wall(0.7,-0.5),Wall(0.7,0.5))
# sde = SDE_VIO(f2,g2,Ws, par)
# axisgrid = (QuinticAxis(Ws[1].pos,Ws[2].pos,71), QuinticAxis(-3.,3.,71))


# @run sdestep = SDEStep(sde, method, 0.1);
# @run PathIntegrationMethod.compute_velocities_to_impact!(step2, sdestep.ID[])

method = Euler()
Δt = 0.1
##

@time PI = PathIntegration(sde, method, Δt, axisgrid..., di_N = 31);
for _ in 1:30
    advance!(PI)
end
# reinit_PI_pdf!(PI)
# advance_till_converged!(PI)

##
begin
    figure(1); clf()
    
    
    res = PI.pdf
    X = [res.axes[1][i] for i in eachindex(res.axes[1]), j in eachindex(res.axes[2])]
    Y = [res.axes[2][j] for i in eachindex(res.axes[1]), j in eachindex(res.axes[2])]
    
    # scatter3D(X, Y, abs.(PI.pdf.p) .|> log10)
    scatter3D(X, Y, PI.pdf.p)
    ax = gca()
    ax.elev = 0
    ax.azim = 80
    # scatter3D(xvs[1],-1 .+ 0.25*xvs[1].^3,zero(xvs[1]))
    # xlim(left=-6,right=6)
    # ylim(bottom=-6,top=6)
    # zlim(bottom=0,top=0.5)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")

    figure(2); clf()
    vs = LinRange_fromaxis(PI.pdf.axes[2],1001)
    plot(vs,PI.(PI.pdf.axes[1][1],vs))
    plot(vs,PI.(PI.pdf.axes[1][2],vs))
    plot(vs,PI.(PI.pdf.axes[1][3],vs))
    plot(vs,PI.(PI.pdf.axes[1][4],vs))
end



##
@run PathIntegration(sde, method, Δt, axisgrid...);


## TO check in debug

IK.x1
step1.x1
step2.steptracer.v_i


idxs = PathIntegrationMethod.dense_idx_it(PI.IK) |> collect;
getindex.(PI.pdf.axes,idxs[61])

sum(abs,stepMX(PI)[61,:])

