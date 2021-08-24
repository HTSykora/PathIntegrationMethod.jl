using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack

function plot_pip(pip; save_as = nothing, top = 0.5)
    figure(1); clf()
    _x = LinRange(xs[1],xs[end],101)
    plot(_x,pip.pdgrid.(_x))
    ylim(bottom=-0.02, top=top)

    if save_as isa String
        savefig(save_as)
    end
end
##

# 1D problem:
f(x,p,t) = x-x^3 + 0.5cos(t)
g(x,p,t) = sqrt(2)

sde = SDE(f,g)
xs = Axis(-5,5,51,interpolation = :chebyshev)

ts = LinRange(0,2π,201)
@time pip = PathIntegrationProblem(sde,ts,xs, precompute=true)#, method = RKMaruyama());

@time for _ in 1:10
    advance!(pip)
end
@time begin
    figure(1); clf()
    _x = LinRange(xs[1],xs[end],101)
    plot(_x,myp.(_x), label="Initial PDF")
    plot(_x, pip.pdgrid.(_x), label="Iteration")
    
    legend()
end

## Create an animation
fname = "./examples/Animfiles/scalar_cubic_"
save_as_i(i) = fname*"$(i).png"

increment = 2
n_frames = 2((length(ts)-1)÷increment)


@time pip = PathIntegrationProblem(sde,ts,xs, precompute=true)#, method = RKMaruyama());
plot_pip(pip,save_as = save_as_i(0))
@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    plot_pip(pip,save_as = save_as_i(i))
end

`ffmpeg -y -framerate 24 -start_number 0 -i $(fname)%d.png -vframes $(n_frames+1) $(fname)anim.mp4` |> run
`rm $(fname)*.png`