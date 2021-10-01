using Revise
using PathIntegrationMethod
using BenchmarkTools
using FFTW

##
function foo!(X, x)
    for i in eachindex(X)
        X[i] = sin(x + 0.01i)
    end
end
function foo!(X, x)
    for i in eachindex(X)
        X[i] = PathIntegrationMethod.normal1D(0,1/10i,x)
    end
end

##
gridaxis = GridAxis(-1,1,11, interpolation = :chebyshev)

integratormethod = QuadGKIntegrator()
integratormethod = TrapezoidalIntegrator()
integratormethod = ClenshawCurtisIntegrator()

temp = zeros(100);
i_N = 10001;
@time di_qgk = DiscreteIntegrator(QuadGKIntegrator(), temp, i_N, gridaxis)
@time di_tri = DiscreteIntegrator(TrapezoidalIntegrator(), temp, i_N, gridaxis)
@time di_cci = DiscreteIntegrator(ClenshawCurtisIntegrator(), temp, i_N, gridaxis)
 
@time di_qgk(foo!)
@time di_tri(foo!)
@time di_cci(foo!)

r_qgk,r_tri,r_cci = getfield.((di_qgk, di_tri, di_cci),:res)
sum(abs,r_qgk .- r_tri)
sum(abs,r_qgk .- r_cci)

@btime di($foo!)


@btime foo!($temp,1.)

function fejer(n)
    # Weights of the Fejer2, Clenshaw-Curtis and Fejer1 quadratures by DFTs
    # n>1. Nodes: x_k = cos(k*pi/n)
    N=1:2:n-1; l=length(N); m=n-l; K=0:m-1;

    n2 = n÷2
    v = zeros(n);
    for k in 0:(n2-1)
        v[k+1] = 2/(1-4k^2)
    end
    v[n2+1] = (n-3)/(2n2 - 1) - 1
    for k in 0:(n2 -1)
        v[n-k] = v[k+1]
    end
    v
    
    g = zeros(n)
    w0cc = 1/(n^2-1+mod(n,2))
    for k in 0:(n2-1)
        g[k+1] = -w0cc
    end
    g[n2+1] = w0cc*((2-mod(n,2))*n-1)
    for k in 0:(n2 -1)
        g[n-k] = g[k+1]
    end
    v,g
    # ifft(v+g)

    v0=[2 ./ N ./ (N .- 2); 1/N[end]; zeros(m);]
    v2=-v0[1:end-1]-v0[end:-1:2]; wf2=ifft(v2);

    g0=-ones(n); g0[l+1]=g0[l+1]+n; g0[m+1]=g0[m+1]+n;
    g=g0./(n^2-1+mod(n,2)); wcc=ifft(v2+g);
end



n = 11
N=1:2:n-1; l=length(N); m=n-l; K=0:m-1;

v0=[2 ./ N ./ (N .- 2); 1/N[end]; zeros(m)]
v2=-v0[1:end-1]-v0[end:-1:2]; wf2=ifft(v2);

g0=-ones(n); g0[l+1]=g0[l+1]+n; g0[m+1]=g0[m+1]+n;
g=g0./(n^2-1+mod(n,2)); wcc=ifft(v2+g);




function mycc(n)
    d = [2/(1-(2k)^2) for k in 0:n÷2+mod(n,2)]
    d[1] = d[1]/2
    d[end] = d[end]/2
    d
end

@time begin
    N1 = 10001; b = 1; a = -1;
    N=N1-1; bma=(b-a);
    c=zeros(2(N1-1));
    Ns = 2:2:N;

    c[1] = 2.;
    for (i,ic) in enumerate(3:2:N1)
        c[ic] = 2/(1-Ns[i]^2)
    end
    for (i,ic) in enumerate(3:2:N1-1)
        c[end-ic+2] = 2/(1-Ns[i]^2)
    end
    (length(f) == length(c)) |> println
    f = ifft(c) .|> real
    w = zeros(N1);
    w[1] = 0.5f[1]
    for i in 2:length(w)-1
        w[i] = f[i]
    end
    w[end] = 0.5f[N1]
    
    w .= w .* bma
    # c[2,2] = 1.;

    xs =  PathIntegrationMethod.chebygrid(N1)
    display(c)
end


# res .= res .* bma

