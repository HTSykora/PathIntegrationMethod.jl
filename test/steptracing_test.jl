using Revise
using PathIntegrationMethod

f1(x,p,t) = x[2] + t
function f2(x,p,t)
    ζ,λ = p
    - x[1] - λ*x[1]^3 - 2*ζ*x[2] + sin(t)
end
function g(x,p,t)
    sqrt(2)
end
x0 = [1.,2.]; t0 = 0.; t1 = 0.1;
par = [0.1,1.]
sde = SDE((f1,f2), g, par)
sdestep = SDEStep(sde,Euler(),x0, similar(x0), t0,t1)



# SDE(DriftTerm(f1,f2),DiffusionTerm(2,2,1,g), par)