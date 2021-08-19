using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack

struct ExcitationFunction{fT,pT}
    f::fT
    p::pT
end
(ef::ExcitationFunction)(p,t) = ef.f(p,t)
(ef::ExcitationFunction)(t) = ef(ef.p,t)
##

function fx(u,p,t)
    f, g = p # f, g, σ
    f(t) + g
end
function gx(u,p,t)
    p[3] # = σ
end
f_per = ExcitationFunction((p,t)-> cos(π*t + p[1]),[0.]) # p = [φ₀]
##

β = π/18; r = 0.3; d = 0.25;
Mg_F = 0.1245/5 * 9.81 # M*g/|F|
g = Mg_F*sin(β) # M*g/|F|*sin(β)
σ = 0.01
p=[f_per, g, σ]

r = 0.3

sde = SDE_Oscillator1D(fx,gx,par = p)

Nᵥ = 31; Nₓ = 31;
v_ax = Axis(-1,1,31)
ts = LinRange(0,2,31)
vi_sde, pdgrid = create_symmetric_VI_PDGrid(sde, d, r, v_ax, Nₓ)

@time pip = PathIntegrationProblem(vi_sde, pdgrid, ts; precompute=true, σ_init = 0.2);

# @run PathIntegrationProblem(vi_sde, pdgrid, ts; precompute=true);


for _ in 1:10
    advance!(pip)
end
begin
    figure(1); clf()
    
    res = pip.pdgrid
    X = [res.xs[1][i] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    Y = [res.xs[2][j] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
   
    scatter3D(X, Y, pip.pdgrid.p)
    # scatter3D(xvs[1],-1 .+ 0.25*xvs[1].^3,zero(xvs[1]))
    # xlim(left=-6,right=6)
    # ylim(bottom=-6,top=6)
    # zlim(bottom=0,top=0.15)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
end



## Create an animation
function export_to_png(pip,i)
    begin
        figure(1); clf()
        res = pip.pdgrid
        X = [res.xs[1][i] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
        Y = [res.xs[2][j] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
       
        scatter3D(X, Y, pip.pdgrid.p)
        # scatter3D(xvs[1],-1 .+ 0.25*xvs[1].^3,zero(xvs[1]))
        # xlim(left=-6,right=6)
        # ylim(bottom=-6,top=6)
        zlim(bottom=0,top=25)
        xlabel(L"x")
        ylabel(L"v")
        zlabel(L"p(x,v)")
        savefig("./examples/Animfiles/1DVIOscillator_$(i).png")
    end
end

increment = 1
n_frames = 5((length(ts)-1)÷increment)
export_to_png(pip,0)
@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    export_to_png(pip,i)
end

`ffmpeg -y -framerate 24 -start_number 0 -i ./examples/Animfiles/1DVIOscillator_%d.png -vframes $(n_frames+1) ./examples/Animfiles/1DVIOscillator_anim.mp4` |> run






# IK = PathIntegrationMethod.computeintegrationmatrix(vi_sde, pdgrid, ts, EulerMaruyama())
# IK.idx₁ .= [1,1]
# @run ik = IK(IK.temp,-1.)
# quadgk!(IK,IK.temp,v_ax[16],v_ax[17])

# PathIntegrationMethod.get_ξ(EulerMaruyama(),vi_sde,ts[2],ts[1],pdgrid.xs[1][2],nothing,nothing,-1)

