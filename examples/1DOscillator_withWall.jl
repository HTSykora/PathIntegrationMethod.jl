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

Nᵥ = 51; Nₓ = 51; # 300 seconds at this resolution
x_lims = (walls[1].pos, length(walls) >1 ? walls[2].pos : 5.)
x_ax = Axis(x_lims...,Nₓ)
v_ax = Axis(-6,6,Nᵥ)
Δt = 0.05

@time pip = PathIntegrationProblem(vi_sde, Δt, x_ax, v_ax; precompute=true, μ_init=[2.,0.], σ_init = [0.1,0.25]);

# @save "./examples/Animfiles/1DOscillator_withWall_initpip.jld2" pip

@time for _ in 1:25
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


##
# Create point-cloud animation
@load "./examples/Animfiles/1DVIOscillator_initpip_noharmonic.jld2"
fname = "./examples/Animfiles/1DVIOscillator_noharmonic_"

increment = 5; step = ceil(Int,6/Δt)
n_frames = (step÷increment)
save_as_i = i-> fname*"$(i).png"
scatter_pip(pip; save_as = save_as_i(0))
@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    scatter_pip(pip; save_as = save_as_i(0))
end

`ffmpeg -y -framerate 24 -start_number 0 -i $(fname)%d.png -vframes $(n_frames+1) $(fname)anim.mp4` |> run

`rm $(fname)*.png` |> run # clean up pngs


# Create surface animation

_xmin, _xmax, _Nx= -0.25,0.25, 101
_vmin, _vmax, _Nv= -1.,1., 101
xrange, vrange =  LinRange(_xmin,_xmax,_Nx), LinRange(_vmin,_vmax,_Nv);
_X = [x for x in xrange, v in vrange]
_V = [v for x in xrange, v in vrange]
_P = similar(_X)

@load "./examples/Animfiles/1DVIOscillator_initpip.jld2"
fname = "./examples/Animfiles/1DVIOscillator_noharmonic_surf_"
save_as_i = (i-> fname*"$(i).png")

increment = 5; step = ceil(Int,12/Δt)
n_frames = (step÷increment)

@time surf_pip!(_P,_X, _V, xrange, vrange, pip)
@time surf_pip!(_P,_X, _V, xrange, vrange, pip, save_as = save_as_i(0))

@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    surf_pip!(_P,_X, _V, xrange, vrange, pip, save_as = save_as_i(i))
end

`ffmpeg -y -framerate 24 -start_number 0 -i $(fname)%d.png -vframes $(n_frames+1) $(fname)anim.mp4` |> run

`rm $(fname)*.png` |> run # clean up pngs

