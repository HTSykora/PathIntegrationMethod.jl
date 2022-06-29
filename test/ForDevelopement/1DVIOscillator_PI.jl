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

W = Wall(x->0.7-0.05x,-1.)
W = Wall(0.7,-0.1)
sde = SDE_VIO(f2,g2,W, par)
axisgrid = (QuinticAxis(W.pos,3.,51), QuinticAxis(-3.05,3.,51))

# Ws = (Wall(0.7,-1.),Wall(0.7,1.))
# sde = SDE_VIO(f2,g2,Ws, par)
# axisgrid = (QuinticAxis(Ws[1].pos,Ws[2].pos,51), QuinticAxis(-3.,3.,51))


# @run sdestep = SDEStep(sde, method, 0.1);
# @run PathIntegrationMethod.compute_velocities_to_impact!(step2, sdestep.ID[])

method = Euler()
Δt = 0.01
##

PI = PathIntegration(sde, method, Δt, axisgrid...);
reinit_PI_pdf!(PI)
for _ in 1:50
    advance!(PI)
end
advance_till_converged!(PI)

##
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

begin
    figure(1); clf()
    
    
    res = PI.pdf
    X = [res.axes[1][i] for i in eachindex(res.axes[1]), j in eachindex(res.axes[2])]
    Y = [res.axes[2][j] for i in eachindex(res.axes[1]), j in eachindex(res.axes[2])]
    
    # scatter3D(X, Y, abs.(PI.pdf.p) .|> log10)
    scatter3D(X, Y, PI.pdf.p)
    ax = gca()
    ax.elev = 20
    ax.azim = -160
    # scatter3D(xvs[1],-1 .+ 0.25*xvs[1].^3,zero(xvs[1]))
    # xlim(left=-6,right=6)
    # ylim(bottom=-6,top=6)
    # zlim(bottom=0,top=0.5)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
end



##
    # @run PathIntegration(sde, method, Δt, axisgrid...);
