using Pkg; Pkg.activate()
using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
using JLD2
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack

##

function fx(u,p,t)
    g = p[1] # f, g, σ
    g
end
function gx(u,p,t)
    p[2] # = σ
end
##

β = π/18; r = 0.3; d = 0.25;
β = π/4; r = 0.3; d = 0.25;
Mg_F = 0.1245/5 * 9.81 # M*g/|F|
g = Mg_F*sin(β) # M*g/|F|*sin(β)
σ = 0.05
p=(g, σ)

r = 0.8

sde = SDE_Oscillator1D(fx,gx,par = p)

Nᵥ = 71; Nₓ = 71; # 300 seconds at this resolution
v_ax = Axis(-1,1,Nᵥ)
Δt = 0.005
vi_sde, pdgrid = create_symmetric_VI_PDGrid(sde, d, r, v_ax, Nₓ,σ_init =  [0.001,0.05])

@time pip = PathIntegrationProblem(vi_sde, pdgrid, Δt; precompute=true);

@save "./examples/Animfiles/1DVIOscillator_initpip_noharmonic.jld2" pip


for _ in 1:500
    advance!(pip)
end
begin
    figure(1); clf()
    
    res = pip.pdgrid
    X = [res.xs[1][i] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    Y = [res.xs[2][j] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    
    scatter3D(X, Y, res.p)
    ax = gca();
    ax.view_init(elev=60, azim = 45)
    # scatter3D(xvs[1],-1 .+ 0.25*xvs[1].^3,zero(xvs[1]))
    # xlim(left=-6,right=6)
    # ylim(bottom=-6,top=6)
    zlim(bottom=0,top=15)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
end


##
# Create an animation
@load "./examples/Animfiles/1DVIOscillator_initpip_noharmonic.jld2"
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
        zlim(bottom=0,top=10)
        xlabel(L"x")
        ylabel(L"v")
        zlabel(L"p(x,v)")
        savefig("./examples/Animfiles/1DVIOscillator_$(i)_noharmonic.png")
    end
end

increment = 5; step = ceil(Int,6/Δt)
n_frames = (step÷increment)
export_to_png(pip,0)
@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    export_to_png(pip,i)
end

`ffmpeg -y -framerate 24 -start_number 0 -i ./examples/Animfiles/1DVIOscillator_%d_noharmonic.png -vframes $(n_frames+1) ./examples/Animfiles/1DVIOscillator_noharmonic_anim.mp4` |> run




_xmin, _xmax, _Nx= -0.25,0.25, 101
_vmin, _vmax, _Nv= -1.,1., 101
xrange, vrange =  LinRange(_xmin,_xmax,_Nx), LinRange(_vmin,_vmax,_Nv);
_X = [x for x in xrange, v in vrange]
_V = [v for x in xrange, v in vrange]
_P = similar(_X)

function export_surf_to_png!(_P, _X, _V,xrange, vrange, pip,i)
    for (j,v) in enumerate(vrange)
        for (i,x) in enumerate(xrange)
            _P[i,j] = pip.pdgrid(x,v)
        end
    end
    figure(1); clf()
    
    plot_surface(_X, _V, _P, cmap=PyPlot.cm.jet)

    ax = gca()
    ax.view_init(elev=60, azim = -90)
    
    zlim(bottom = 0, top = 20)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
    savefig("./examples/Animfiles/1DVIOscillator_noharmonic_surf_$(i).png")
end

@load "./examples/Animfiles/1DVIOscillator_initpip.jld2"
increment = 5; step = ceil(Int,12/Δt)
n_frames = (step÷increment)
@time export_surf_to_png!(_P,_X, _V, xrange, vrange, pip, 0)

@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    export_surf_to_png!(_P,_X, _V, xrange, vrange, pip, i)
end

`ffmpeg -y -framerate 24 -start_number 0 -i ./examples/Animfiles/1DVIOscillator_noharmonic_surf_%d.png -vframes $(n_frames+1) ./examples/Animfiles/1DVIOscillator_noharmonic_surf_anim.mp4` |> run

##
# IK = PathIntegrationMethod.computeintegrationmatrix(vi_sde, pdgrid, ts, EulerMaruyama())


# PathIntegrationMethod.update_idx1!(IK,(3,1))
# PathIntegrationMethod.get_integ_limits(IK)
# @run ik = IK(IK.temp,-1.)
# quadgk!(IK,IK.temp,v_ax[16],v_ax[17])

# PathIntegrationMethod.get_ξ(EulerMaruyama(),vi_sde,ts[2],ts[1],pdgrid.xs[1][2],nothing,nothing,-1)

