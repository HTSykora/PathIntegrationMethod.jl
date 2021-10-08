using Pkg; Pkg.activate()
using Revise
using PathIntegrationMethod
using BenchmarkTools
##

f1(x,p,t) = x[2]^2 + x[1]^2 + x[3]^2 + t
function f2(x,p,t)
    ζ,λ = p
    - x[1] - λ*x[1]^3 - 2*ζ*x[2] + sin(t) + λ*x[3]^3 + 2ζ*x[3]
end
function f3(x,p,t)
    ζ,λ = p
    - x[1] - λ*x[1]^3 - 2*ζ*x[2] + cos(t) - λ*x[3]^2
end
function g(x,p,t)
    sqrt(2)
end

method = Euler()
method = RK4()
##
x0 = [1., 2., 3.]; t0 = 0.6; t1 = t0 + 0.1;
par = [0.1, 1.1]
sde = SDE((f1, f2, f3), g, par)

@time sdestep = SDEStep(sde, method, x0, similar(x0), t0, t1);

@time PathIntegrationMethod.eval_driftstep!(sdestep)
x0_ref = deepcopy(sdestep.x0)
x1_ref = deepcopy(sdestep.x1)

sdestep.x0[1] = x1_ref[1]
sdestep.x0[2] = x1_ref[2]
sdestep.x1[3] = sdestep.x0[3]

##
@time PathIntegrationMethod.compute_missing_states_driftstep!(sdestep)
# @btime PathIntegrationMethod.compute_missing_states_driftstep!($sdestep)
reduce(&, x0_ref .≈ sdestep.x0)
reduce(&, x1_ref .≈ sdestep.x1)


@run PathIntegrationMethod.compute_missing_states_driftstep!(sdestep)


# RuntimeGeneratedFunctions.RuntimeGeneratedFunction{(:ˍ₋out, :x, :y, :par, :t0, :t1), Symbolics.var"#_RGF_ModTag", Symbolics.var"#_RGF_ModTag", (0xdfc6b714, 0x6f54baa1, 0x981b40b7, 0x838bcf85, 0x4349bad8)}(quote
#     #= /home/hts2000/.julia/packages/SymbolicUtils/vKBkE/src/code.jl:284 =#
#     #= /home/hts2000/.julia/packages/SymbolicUtils/vKBkE/src/code.jl:285 =#
#     begin
#         #= /home/hts2000/.julia/packages/Symbolics/mFWWM/src/build_function.jl:373 =#
#         #= /home/hts2000/.julia/packages/SymbolicUtils/vKBkE/src/code.jl:331 =# @inbounds begin
#                 #= /home/hts2000/.julia/packages/SymbolicUtils/vKBkE/src/code.jl:327 =#
#                 ˍ₋out[1] = (+)((/)((+)((getindex)(x, 1), (*)((+)(t1, (*)(-1, t0)), (+)((^)((getindex)(x, 1), 2), (^)((getindex)(x, 3), 2), t0, (^)((getindex)(x, 2), 2))), (*)(-1, (getindex)(y, 1)), (/)((*)(2, (+)((getindex)(y, 2), (/)((*)(-1, (+)((*)(-1, t1), t0), (+)(-1, (*)(-3, (getindex)(par, 2), (^)((getindex)(x, 1), 2))), (+)((*)(-1, (getindex)(x, 1)), (*)(-1, (+)(t1, (*)(-1, t0)), (+)((^)((getindex)(x, 1), 2), (^)((getindex)(x, 3), 2), t0, (^)((getindex)(x, 2), 2))), (getindex)(y, 1))), (+)(-1, (*)(2, (+)((*)(-1, t1), t0), (getindex)(x, 1)))), (*)(-1, (getindex)(x, 2)), (*)(-1, (+)((+)((sin)(t0), (*)((getindex)(par, 2), (^)((getindex)(x, 3), 3)), (*)(-1, (getindex)(x, 1)), (*)(-2, (getindex)(x, 2), (getindex)(par, 1))), (+)((*)(-1, (getindex)(par, 2), (^)((getindex)(x, 1), 3)), (*)(2, (getindex)(x, 3), (getindex)(par, 1)))), (+)(t1, (*)(-1, t0)))), (+)((*)(-1, t1), t0), (getindex)(x, 2)), (+)(-1, (/)((*)(-2, (^)((+)((*)(-1, t1), t0), 2), (getindex)(x, 2), (+)(-1, (*)(-3, (getindex)(par, 2), (^)((getindex)(x, 1), 2)))), (+)(-1, (*)(2, (+)((*)(-1, t1), t0), (getindex)(x, 1)))), (*)(-2, (+)((*)(-1, t1), t0), (getindex)(par, 1))))), (+)(-1, (*)(2, (+)((*)(-1, t1), t0), (getindex)(x, 1)))), (getindex)(x, 1))
#                 ˍ₋out[2] = (+)((/)((+)((*)(-1, (getindex)(y, 2)), (/)((*)((+)((*)(-1, t1), t0), (+)(-1, (*)(-3, (getindex)(par, 2), (^)((getindex)(x, 1), 2))), (+)((*)(-1, (getindex)(x, 1)), (*)(-1, (+)(t1, (*)(-1, t0)), (+)((^)((getindex)(x, 1), 2), (^)((getindex)(x, 3), 2), t0, (^)((getindex)(x, 2), 2))), (getindex)(y, 1))), (+)(-1, (*)(2, (+)((*)(-1, t1), t0), (getindex)(x, 1)))), (getindex)(x, 2), (*)((+)((+)((sin)(t0), (*)((getindex)(par, 2), (^)((getindex)(x, 3), 3)), (*)(-1, (getindex)(x, 1)), (*)(-2, (getindex)(x, 2), (getindex)(par, 1))), (+)((*)(-1, (getindex)(par, 2), (^)((getindex)(x, 1), 3)), (*)(2, (getindex)(x, 3), (getindex)(par, 1)))), (+)(t1, (*)(-1, t0)))), (+)(-1, (/)((*)(-2, (^)((+)((*)(-1, t1), t0), 2), (getindex)(x, 2), (+)(-1, (*)(-3, (getindex)(par, 2), (^)((getindex)(x, 1), 2)))), (+)(-1, (*)(2, (+)((*)(-1, t1), t0), (getindex)(x, 1)))), (*)(-2, (+)((*)(-1, t1), t0), (getindex)(par, 1)))), (getindex)(x, 2))
#                 #= /home/hts2000/.julia/packages/SymbolicUtils/vKBkE/src/code.jl:329 =#
#                 nothing
#             end
#     end
# end)