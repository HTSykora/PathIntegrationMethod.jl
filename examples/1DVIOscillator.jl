include("./example_init.jl")

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
σ = 0.1
p=(f_per, -g, σ)

r = 0.3

sde = SDE_Oscillator1D(fx,gx,par = p)

Nᵥ = 51; Nₓ = 51; # 300 seconds at this resolution
v_ax = Axis(-1.5,1.5,Nᵥ)
ts = LinRange(0,2,51)
vi_sde, pdgrid = create_symmetric_VI_PDGrid(sde, d, r, v_ax, Nₓ,σ_init = [0.025,0.25])

@time pip = PathIntegrationProblem(vi_sde, pdgrid, ts; precompute=true);

@save "./examples/Animfiles/1DVIOscillator_initpip.jld2" pip
# @run PathIntegrationProblem(vi_sde, pdgrid, ts; precompute=true);


for _ in 1:2
    advance!(pip)
end
scatter_pip(pip; elev=60, azim = 45., top = 100)

##
# Create an animation
@load "./examples/Animfiles/1DVIOscillator_initpip.jld2"
fname = "./examples/Animfiles/1DVIOscillator_"
save_as_i(i) = fname*"$(i).png"


increment = 1
n_frames = 5((length(ts)-1)÷increment)
scatter_pip(pip; elev=60, azim = 45., top = 100, save_as = save_as_i(0))
@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    scatter_pip(pip; elev=60, azim = 45., top = 100, save_as = save_as_i(i))
end

`ffmpeg -y -framerate 24 -start_number 0 -i $(fname)%d.png -vframes $(n_frames+1) $(fname)_anim.mp4` |> run
`rm $(fname)*.png` |> run # clean up pngs




_xmin, _xmax, _Nx= -0.25,0.25, 101
_vmin, _vmax, _Nv= -1.5,1.5, 101
xrange, vrange =  LinRange(_xmin,_xmax,_Nx), LinRange(_vmin,_vmax,_Nv);
_X = [x for x in xrange, v in vrange]
_V = [v for x in xrange, v in vrange]
_P = similar(_X)

@load "./examples/Animfiles/1DVIOscillator_initpip.jld2"
fname = "./examples/Animfiles/1DVIOscillator_surf_"
save_as_i(i) = fname*"$(i).png"

increment = 1
n_frames = 5((length(ts)-1)÷increment)
@time surf_pip!(_P,_X, _V, xrange, vrange, pip, save_as = save_as_i(0))

@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    surf_pip!(_P,_X, _V, xrange, vrange, pip, save_as = save_as_i(i))
end
`ffmpeg -y -framerate 24 -start_number 0 -i $(fname)%d.png -vframes $(n_frames+1) $(fname)_anim.mp4` |> run
`rm $(fname)*.png` |> run # clean up pngs
