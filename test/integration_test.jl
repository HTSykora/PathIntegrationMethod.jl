using PathIntegrationMethod

##
function foo!(X, x)
    for i in eachindex(X)
        X[i] = sin(x + 0.01i)
    end
end
function foo!(X, x)
    for i in eachindex(X)
        X[i] = PathIntegrationMethod.normal1D(i/10length(X),1/i,x)
    end
end

##
gridaxis = ChebyshevAxis(-3,3,7)

integratormethod = QuadGKIntegrator()
integratormethod = TrapezoidalIntegrator()
integratormethod = NewtonCotesIntegrator(2)
integratormethod = NewtonCotesIntegrator(3)
integratormethod = ClenshawCurtisIntegrator()
integratormethod = GaussLegendreIntegrator()
integratormethod = GaussRadauIntegrator()
integratormethod = GaussLobattoIntegrator()

temp = zeros(100);
i_N = 2001;
# @time 
di_qgk = DiscreteIntegrator(QuadGKIntegrator(), temp, i_N, gridaxis)
# @time 
di_tri = DiscreteIntegrator(TrapezoidalIntegrator(), temp, i_N, gridaxis)
# @time 
di_nc2 = DiscreteIntegrator(NewtonCotesIntegrator(2), temp, i_N, gridaxis)
# @time 
di_nc3 = DiscreteIntegrator(NewtonCotesIntegrator(3), temp, i_N, gridaxis)
# @time 
di_cci = DiscreteIntegrator(ClenshawCurtisIntegrator(), temp, i_N, gridaxis)
# @time 
di_glei = DiscreteIntegrator(GaussLegendreIntegrator(), temp, i_N, gridaxis)
# @time 
di_gri = DiscreteIntegrator(GaussRadauIntegrator(), temp, i_N, gridaxis)
# @time 
di_gloi = DiscreteIntegrator(GaussLobattoIntegrator(), temp, i_N, gridaxis)
 
# @time 
di_qgk(foo!)
# @time 
di_tri(foo!)
# @time 
di_nc2(foo!)
# @time 
di_nc3(foo!)
# @time 
di_cci(foo!)
# @time 
di_glei(foo!)
# @time 
di_gri(foo!)
# @time 
di_gloi(foo!)

r_qgk,r_tri, r_nc2, r_nc3, r_cci, r_glei, r_gri, r_gloi = getfield.((di_qgk, di_tri, di_nc2, di_nc3, di_cci, di_glei, di_gri, di_gloi),:res)
testresults = BitArray(undef,0)
push!(testresults, sum(abs,r_qgk .- r_tri) < 1e-7)
push!(testresults, sum(abs,r_qgk .- r_nc2) < 1e-12)
push!(testresults, sum(abs,r_qgk .- r_nc3) < 1e-6)
push!(testresults, sum(abs,r_qgk .- r_cci) < 1e-12)
push!(testresults, sum(abs,r_qgk .- r_glei) < 1e-12)
push!(testresults, sum(abs,r_qgk .- r_gri) < 1e-12)
push!(testresults, sum(abs,r_qgk .- r_gloi) < 1e-12)
reduce(&,testresults)