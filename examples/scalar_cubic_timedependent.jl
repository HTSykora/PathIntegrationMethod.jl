using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack
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
function export_to_png(pip,i)
    begin
        figure(1); clf()
        _x = LinRange(xs[1],xs[end],101)
        plot(_x,pip.pdgrid.(_x), label="Iteration")
        
        ylim(bottom=-0.02,top=0.5)
        savefig("./examples/Animfiles/scalar_cubic_$(i).png")
    end
end


increment = 2
n_frames = 2((length(ts)-1)÷increment)
@time pip = PathIntegrationProblem(sde,ts,xs, precompute=true)#, method = RKMaruyama());
export_to_png(pip,0)
@time for i in 1:n_frames
    for _ in 1:increment
        advance!(pip)
    end
    export_to_png(pip,0)
end

`ffmpeg -y -framerate 24 -start_number 0 -i ./examples/Animfiles/scalar_cubic_%d.png -vframes $(n_frames+1) ./examples/Animfiles/scala_cubic_anim.mp4` |> run