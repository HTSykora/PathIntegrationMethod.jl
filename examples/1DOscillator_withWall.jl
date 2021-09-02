include("./example_init.jl")
##

function fx(u,p,t)
    ζ = p[1] # f, g, σ
    -2ζ*u[2] - u[1]
end
function gx(u,p,t)
    p[2] # = σ
end
##

ζ = 0.1; σ = 1.
p=[ζ, σ]

r = 0.7
# walls = (Wall(r,0.0),Wall(r,3.))
walls = (Wall(r,0.),)
sde = SDE_Oscillator1D(fx,gx,par = p)
vi_sde = SDE_VI_Oscillator1D(sde,walls)

Nᵥ = 71; Nₓ = 71; 
x_lims = (walls[1].pos, length(walls) >1 ? walls[2].pos : 6.)
x_ax = Axis(x_lims...,Nₓ)
v_ax = Axis(-6,6,Nᵥ)
Δt = 0.05

@time pip = PathIntegrationProblem(vi_sde, Δt, x_ax, v_ax; precompute = true, μ_init = [2.,0.], σ_init = [0.1,0.25]);

@save "./examples/Animfiles/1DOscillator_withWall_initpip.jld2" pip

@time for _ in 1:100
    advance!(pip)
end
scatter_pip(pip; elev=45, azim = 145, top = 0.8)
# scatter3D(fill(pip.pdgrid.xs[1][1],Nᵥ),pip.pdgrid.xs[2],pip.pdgrid[1,:])
# scatter3D(fill(pip.pdgrid.xs[1][2],Nᵥ),pip.pdgrid.xs[2],pip.pdgrid[2,:])
# scatter3D(fill(pip.pdgrid.xs[1][3],Nᵥ),pip.pdgrid.xs[2],pip.pdgrid[3,:])
begin
    figure(2); clf();
    for i in 1:10
        plot(pip.pdgrid.xs[2],pip.pdgrid[i,:], label=i)
    end
    legend()
end

## Distribution of each state variable
@time for _ in Δt:Δt:(6+sqrt(eps())) # == 120Δt
    advance!(pip)
end
x_dist = sum(pip.pdgrid[:,i]*wt for (i,wt) in enumerate(pip.pdgrid.xs[2].wts))
v_dist = sum(pip.pdgrid[i,:]*wt for (i,wt) in enumerate(pip.pdgrid.xs[1].wts))
p_x=PDGrid(x_dist,pip.pdgrid.xs[1])
p_v=PDGrid(v_dist,pip.pdgrid.xs[2])

begin
    fig = figure("1DOscillator_withWall_Distributions", figsize=(10/2.54,7/2.54)); clf()
    xvals = LinRange(p_x.xs[1][1], p_x.xs[1][end], 501)
    vvals = LinRange(p_v.xs[1][1], p_v.xs[1][end], 501)

    subplot(211)
    ax = gca()
    ax.plot(xvals,p_x.(xvals))
    xlabel(L"x", labelpad=1, fontsize=10)
    ylabel(L"p_x(x)", labelpad=3, fontsize=10)
    ylim(top=1.1)
    subplot(212)
    ax = gca()
    ax.plot(vvals,p_v.(vvals))
    xlabel(L"v", labelpad=1, fontsize=10)
    ylabel(L"p_v(v)", labelpad=3, fontsize=10)
    ylim(top=0.45)

    tight_layout(pad=0.2,rect=[0.005, 0.001, 0.999, 0.999])
    fig.canvas.draw() # Update the figure
end
##
# Create point-cloud animation
@load "./examples/Animfiles/1DOscillator_withWall_initpip.jld2"
fname = "./examples/Animfiles/1DOscillator_withWall_"

increment = 2; step = ceil(Int,24/Δt)
n_frames = (step÷increment)
save_as_i = i-> fname*"$(i).png"
scatter_pip(pip; save_as = save_as_i(0), top = 0.7)
@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    scatter_pip(pip; save_as = save_as_i(i), top = 0.7)
end

begin
    `ffmpeg -y -framerate 24 -start_number 0 -i $(fname)%d.png -vframes $(n_frames+1) $(fname)anim.mp4` |> run
    [`rm $(fname)$i.png` |> run for i in 0:n_frames] # clean up pngs
end
# Create surface animation

_xmin, _xmax, _Nx= x_lims..., 501
_vmin, _vmax, _Nv= v_ax[1], v_ax[end], 501
xrange, vrange =  LinRange(_xmin,_xmax,_Nx), LinRange(_vmin,_vmax,_Nv);
_X = [x for x in xrange, v in vrange]
_V = [v for x in xrange, v in vrange]
_P = similar(_X)

@load "./examples/Animfiles/1DOscillator_withWall_initpip.jld2"
fname = "./examples/Animfiles/1DVIOscillator_noharmonic_surf_"
save_as_i = (i-> fname*"$(i).png")

increment = 2; step = ceil(Int,12/Δt)
n_frames = (step÷increment)

@time surf_pip!(_P,_X, _V, xrange, vrange, pip)
@time surf_pip!(_P,_X, _V, xrange, vrange, pip, save_as = save_as_i(0), top = 0.7)

@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    surf_pip!(_P,_X, _V, xrange, vrange, pip, save_as = save_as_i(i), top = 0.7)
end
begin
    `ffmpeg -y -framerate 24 -start_number 0 -i $(fname)%d.png -vframes $(n_frames+1) $(fname)anim.mp4` |> run
    [`rm $(fname)$i.png` |> run for i in 0:n_frames] # clean up pngs
end

