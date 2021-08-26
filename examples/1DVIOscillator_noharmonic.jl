include("./example_init.jl")
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
Mg_F = 0.1245/5 * 9.81 # M*g/|F|
g = Mg_F*sin(β) # M*g/|F|*sin(β)
σ = 0.05
p=(g, σ)

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
scatter_pip(pip; elev=60, azim = 45)


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

@time surf_pip!(_P,_X, _V, xrange, vrange, pip, save_as = save_as_i(0))

@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    surf_pip!(_P,_X, _V, xrange, vrange, pip, save_as = save_as_i(i))
end

`ffmpeg -y -framerate 24 -start_number 0 -i $(fname)%d.png -vframes $(n_frames+1) $(fname)anim.mp4` |> run

`rm $(fname)*.png` |> run # clean up pngs
