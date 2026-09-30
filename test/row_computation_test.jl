module RowComputationTest
using PathIntegrationMethod, Test
const PIM = PathIntegrationMethod

# Scalar system: dx = (x - x³) dt + √2 dW
f(x,p,t) = x[1] - p[1]*x[1]^3
g(x,p,t) = sqrt(2)
# Duffing oscillator: dx = v dt, dv = (-2ζv + x - λx³) dt + σ dW
f1(x,p,t) = x[2]
f2(x,p,t) = -2p[1]*x[2] + x[1] - p[2]*x[1]^3
g2(x,p,t) = p[3]
# Scalar system with periodic forcing: dx = (x - x³ + A sin(2πt/T)) dt + √2 dW
fp(x,p,t) = x[1] - x[1]^3 + p[1]*sin(2π*t/p[2])

# (label, PathIntegration constructor with parameters par and keyword arguments, parameters)
systems = (("d=1 cubic", (par; kw...) -> PathIntegration(SDE(f, g, par), RK2(), 0.05, CubicAxis(-3., 3., 41); kw...), [1.0]),
           ("d=1 Chebyshev", (par; kw...) -> PathIntegration(SDE(f, g, par), RK2(), 0.05, ChebyshevAxis(-3., 3., 31); kw...), [1.0]),
           ("d=2 quintic", (par; kw...) -> PathIntegration(SDE((f1, f2), g2, par), Euler(), 0.05, QuinticAxis(-4., 4., 15), QuinticAxis(-4., 4., 15); kw...), [0.5, 0.25, 1.0]),
           ("d=2 quintic, dense step matrix", (par; kw...) -> PathIntegration(SDE((f1, f2), g2, par), RK4(), 0.05, QuinticAxis(-4., 4., 15), QuinticAxis(-4., 4., 15); stepMXtype = DenseMX(), kw...), [0.5, 0.25, 1.0]),
           ("d=1 QuadGK", (par; kw...) -> PathIntegration(SDE(f, g, par), RK4(), 0.05, CubicAxis(-3., 3., 21); discreteintegrator = QuadGKIntegrator(), kw...), [1.0]),
           ("time-periodic", (par; kw...) -> PathIntegration(SDE(fp, g, par), Euler(), collect(range(0, 0.3, length = 4)), CubicAxis(-3., 3., 41); kw...), [1.0, 0.5]))

# more chunks than rows: some chunks are empty
rowcomputations = (ThreadedRowComputation(1), ThreadedRowComputation(3), ThreadedRowComputation(), ThreadedRowComputation(50),
                   BatchRowComputation(3), BatchRowComputation(), BatchRowComputation(50))

stepmatrices(PI) = [Matrix(S) for S in PI.stepMX]
# the stored matrices, including the sparsity structure
same_storage(A, B) = A == B
same_storage(A::PIM.SparseMatrixCSC, B::PIM.SparseMatrixCSC) = A.colptr == B.colptr && A.rowval == B.rowval && A.nzval == B.nzval
same_stepMX(PI1, PI2) = all(same_storage(PIM.storage_matrix(S1), PIM.storage_matrix(S2)) for (S1, S2) in zip(PI1.stepMX, PI2.stepMX))

@testset "Parallel row computation equals the serial one: $lbl" for (lbl, new_PI, par) in systems
    PI_serial = new_PI(copy(par); rowcomputation = SerialRowComputation())
    @test length(PI_serial.stepMX) == (lbl == "time-periodic" ? 3 : 1)
    for rowcomputation in rowcomputations
        PI = new_PI(copy(par); rowcomputation)
        @test same_stepMX(PI, PI_serial)
        @test PI.stepMX_wts == PI_serial.stepMX_wts
    end
end

@testset "recompute_PI! with parallel row computation: $lbl" for (lbl, new_PI, par) in systems[[1, 3, 6]]
    par_new = 0.8 .* par
    PI_new = new_PI(copy(par_new); rowcomputation = SerialRowComputation())
    for rowcomputation in (ThreadedRowComputation(3), BatchRowComputation(3))
        PI = new_PI(copy(par); rowcomputation)
        advance!(PI)
        recompute_PI!(PI; par = par_new, Q_reinit_pdf = true) # the copies of the integration kernel use the new parameters
        @test same_stepMX(PI, PI_new)
        @test PI.stepMX_wts == PI_new.stepMX_wts
        advance!(PI); advance!(PI_new)
        @test PI.pdf.p == PI_new.pdf.p
        reinit_PI_pdf!(PI_new)
    end
    # the row computation can be changed at recomputation
    PI = new_PI(copy(par_new); rowcomputation = SerialRowComputation())
    recompute_stepMX!(PI; rowcomputation = ThreadedRowComputation(4))
    @test same_stepMX(PI, PI_new)
end

@testset "Row computation types" begin
    for RC in (ThreadedRowComputation, BatchRowComputation)
        @test RC().N_threads == Threads.nthreads()
        @test RC(4).N_threads == 4
        @test RC(N_threads = 5).N_threads == 5
        @test_throws ArgumentError RC(0)
    end
    @test PIM.default_rowcomputation() == (Threads.nthreads() > 1 ? ThreadedRowComputation() : SerialRowComputation())
    PI = systems[1][2](copy(systems[1][3]))
    @test PI.IK.kwargs.rowcomputation == PIM.default_rowcomputation()
    @test PIM.chunk_ranges(10, 3) == [1:3, 4:6, 7:10]
    @test PIM.chunk_ranges(2, 3) == [1:0, 1:1, 2:2]
end

@testset "Nested in a threaded loop" begin
    # `@threads :static` cannot be nested: the chunks are computed one after the other
    lbl, new_PI, par = systems[3]
    S_serial = stepmatrices(new_PI(copy(par); rowcomputation = SerialRowComputation()))
    PIs = [new_PI(copy(par); rowcomputation = SerialRowComputation()) for _ in 1:2]
    for rowcomputation in (ThreadedRowComputation(3), BatchRowComputation(3))
        Threads.@threads :static for PI in PIs
            recompute_stepMX!(PI; rowcomputation)
        end
        @test all(stepmatrices(PI) == S_serial for PI in PIs)
    end
end

end
