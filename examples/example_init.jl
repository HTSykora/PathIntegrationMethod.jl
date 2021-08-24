using Pkg; Pkg.activate()
using Revise, BenchmarkTools
using PathIntegrationMethod
using PyPlot, LaTeXStrings; pygui(true);
using JLD2
PyPlot.rc("text", usetex=true);
py_colors=PyPlot.PyDict(PyPlot.matplotlib."rcParams")["axes.prop_cycle"].by_key()["color"];

using QuadGK, Arpack


function scatter_pip(pip; azim = 45, elev = 60, top = 10, save_as = nothing, fig_id = 1)
    figure(fig_id); clf()
    res = pip.pdgrid
    X = [res.xs[1][i] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    Y = [res.xs[2][j] for i in eachindex(res.xs[1]), j in eachindex(res.xs[2])]
    
    scatter3D(X, Y, pip.pdgrid.p)
    ax = gca()
    ax.view_init(elev = elev, azim = azim)
    
    zlim(bottom=0,top=10)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
    if save_as isa String
        savefig(save_as)
    end
end
function surf_pip!(_P, _X, _V,xrange, vrange, pip; azim = 45, elev = 60, top = 20, save_as = nothing, fig_id = 1)
    for (j,v) in enumerate(vrange)
        for (i,x) in enumerate(xrange)
            _P[i,j] = pip.pdgrid(x,v)
        end
    end
    figure(1); clf()
    
    plot_surface(_X, _V, _P, cmap=PyPlot.cm.jet)

    ax = gca()
    ax.view_init(elev = elev, azim = azim)
    
    zlim(bottom = 0, top = top)
    xlabel(L"x")
    ylabel(L"v")
    zlabel(L"p(x,v)")
    if save_as isa String
        savefig(save_as)
    end
end
##