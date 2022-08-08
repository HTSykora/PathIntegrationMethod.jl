using PathIntegrationMethod


f1(x,p,t) = x[2]^2 + x[1]^2 + x[3]^2 + t
function f2(x,p,t)
    ζ,λ = p
    - x[1] - λ*x[1]^3 - 2*ζ*x[2] + sin(t) + λ*x[3]^3 + 2ζ*x[3]
end
function f3(x,p,t)
    ζ,λ = p
    - x[1] - λ*x[1]^3 - 2*ζ*x[2] + cos(t) - λ*x[3]^2
end
function g1(x,p,t)
    sqrt(2)
end
function g2(x,p,t)
    x[1] + sqrt(2)
end

##
x0 = [1.,2.,3.]; t0 = 0.6; t1 = t0 + 0.1;
par = [0.1,1.1]

# Single noise source
sde = SDE(f1, g1, par)
sde = SDE((f1,), (g1,), par)
# sde = SDE(f1, (g1,g2), par)
sde = SDE((f1,f2,f3), g1, par)
sde = SDE((f1,f2,f3), (g1,g2), par)

true