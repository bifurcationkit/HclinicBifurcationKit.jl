# using Revise
# using Plots
using Test
using BifurcationKit, LinearAlgebra, Setfield, ForwardDiff, Parameters, HclinicBifurcationKit
const BK = BifurcationKit
norminf(x) = norm(x, Inf)

####################################################################################################
# test for the Bogdanov-Takens normal form
function Fbt!(dx, x, p, t=0)
	dx[1] = x[2]
	dx[2]= p.β1 + p.β2 * x[2] + p.a * x[1]^2 + p.b * x[1] * x[2]
	dx
end
Fbt(x,p) = Fbt!(similar(x),x,p)
par = (β1 = -0.01, β2 = -0.1, a = 1., b = -1.)
prob  = BK.BifurcationProblem(Fbt, [0.01, 0.01], par, (@lens _.β1))
opt_newton = NewtonPar(tol = 1e-9, maxIter = 40, verbose = false)
opts_br = ContinuationPar(dsmin = 0.001, dsmax = 0.01, ds = 0.01, pMax = 0.5, pMin = -0.5, detectBifurcation = 3, nev = 2, newtonOptions = opt_newton, maxSteps = 100, nInversion = 8, tolBisectionEigenvalue = 1e-8, dsminBisection = 1e-9, saveSolEveryStep = 1)

br = continuation(prob, PALC(), opts_br; bothside = true, verbosity = 0)

sn_codim2 = continuation(br, 2, (@lens _.β2), ContinuationPar(opts_br, detectBifurcation = 1, saveSolEveryStep = 1, maxSteps = 40) ;
	detectCodim2Bifurcation = 2,
	updateMinAugEveryStep = 1,
	)
@test sn_codim2.specialpoint[1].type == :bt
@test sn_codim2.specialpoint[1].param ≈ 0 atol = 1e-6
@test length(unique(sn_codim2.BT)) == length(sn_codim2)

hopf_codim2 = continuation(br, 3, (@lens _.β2), ContinuationPar(opts_br, detectBifurcation = 1, saveSolEveryStep = 1, maxSteps = 40, maxBisectionSteps = 25) ; plot = false, verbosity = 0,
	detectCodim2Bifurcation = 2,
	updateMinAugEveryStep = 1,
	bothside = true,
	)

@test length(hopf_codim2.specialpoint) == 3
@test hopf_codim2.specialpoint[2].type == :bt
@test hopf_codim2.specialpoint[2].param ≈ 0 atol = 1e-6
@test length(unique(hopf_codim2.BT)) == length(hopf_codim2)-1
# plot(sn_codim2, hopf_codim2, branchlabel = ["Fold", "Hopf"])

# refine BT point
btpt = newton(sn_codim2, 1; startWithEigen = true)

btpt = getNormalForm(sn_codim2, 1; nev = 2, autodiff = false)
show(btpt)
BK.type(btpt)
@test norm(btpt.nf.b * sign(sum(btpt.ζ[1])) - par.b, Inf) < 1e-5
@test norm(btpt.nf.a * sign(sum(btpt.ζ[1])) - par.a, Inf) < 1e-5
@test isapprox(abs.(btpt.ζ[1]), [1, 0])
@test isapprox(abs.(btpt.ζ[2]), [0, 1];rtol = 1e-6)
@test isapprox(abs.(btpt.ζ★[1]), [1, 0];rtol = 1e-6)

@test isapprox(btpt.nfsupp.K2, [0, 0]; atol = 1e-5)
@test isapprox(btpt.nfsupp.d, 0; atol = 1e-3)
@test isapprox(btpt.nfsupp.e, 0; atol = 1e-3)
@test isapprox(btpt.nfsupp.a1, 0; atol = 1e-3)
@test isapprox(btpt.nfsupp.b1, 0; atol = 1e-3)

btpt1 = getNormalForm(sn_codim2, 1; nev = 2, autodiff = false)
@test mapreduce(isapprox, &, btpt.nf, btpt1.nf)
@test mapreduce(isapprox, &, btpt.nfsupp, btpt1.nfsupp)

HC = BK.predictor(btpt, Val(:HopfCurve), 0.)
	HC.hopf(0.)
SN = BK.predictor(btpt, Val(:FoldCurve), 0.)
Hom = BK.predictor(btpt, Val(:HomoclinicCurve), 0.)
	Hom.orbit(0,0)

# plot(sn_codim2, branchlabel = ["Fold"], vars = (:β1, :β2))
# 	_S = LinRange(-0.06, 0.06, 1000)
# 	plot!([HC.hopf(s)[1] for s in _S], [HC.hopf(s)[2] for s in _S], linewidth=5, label = "Hpred")
# 	plot!([SN.fold(s)[1] for s in _S], [SN.fold(s)[2] for s in _S], linewidth=5, label = "SNpred")
# 	_S = LinRange(-0.25, 0.25, 1000)
# 	plot!([Hom.α(s)[1] for s in _S], [Hom.α(s)[2] for s in _S], linewidth=5, label = "Hom")
#
# 	plot!(hopf_codim2, branchlabel = ["Hopf"], vars = (:β1, :β2), color = :black)
# 	xlims!(-0.05, 0.001)


# plot of the homoclinic orbit
hom1 = [Hom.orbit(t,0.1)[1] for t in LinRange(-1000, 1000, 10000)]
hom2 = [Hom.orbit(t,0.1)[2] for t in LinRange(-1000, 1000, 10000)]
# plot(hom1, hom2)

########################
optn_hom = NewtonPar(verbose = false, tol = 1e-10)
optc_hom = ContinuationPar(newtonOptions = optn_hom, ds = -0.0001, dsmin = 1e-5, plotEveryStep = 1,maxSteps = 100, detectBifurcation = 0, saveSolEveryStep = 1)

########################
using OrdinaryDiffEq

function plotHom(x,p;k...)
	𝐇𝐨𝐦 = p.prob
	par0 = set(BK.getParams(𝐇𝐨𝐦), BK.getLens(𝐇𝐨𝐦), x.x[end][1])
	par0 = set(par0, p.lens, p.p)
	sol = getHomoclinicOrbit(𝐇𝐨𝐦, x, par0)
	m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
	plot!(sol.t, sol[:,:]',subplot=3, markersize = 1, marker=m)
end

br_hom = continuation(prob, btpt,
	PeriodicOrbitOCollProblem(20, 4; meshadapt = true),
	PALC(tangent = Bordered()),
	setproperties(optc_hom, maxSteps = 200, saveSolEveryStep = 1, dsmax = 3e-1, plotEveryStep = 3, detectEvent=2);
	amplitude = 5e-3, ϵ0 = 1e-3,
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)), # WORK BEST
	# freeparams = ((@lens _.T),),
	verbosity = 0, plot = false,
	callbackN = BK.cbMaxNorm(1e1),
	plotSolution = plotHom,
	normC = norminf
	)

_sol = getHomoclinicOrbit(br_hom.prob.VF.F, br_hom.sol[end].x, BK.setParam(br_hom,  br_hom.sol[end].p))

# _S = LinRange(-0.06, 0.06, 1000)
# 	plot([HC.hopf(s)[1] for s in _S], [HC.hopf(s)[2] for s in _S], linewidth=5, label = "Hpred")
# 	plot!([SN.fold(s)[1] for s in _S], [SN.fold(s)[2] for s in _S], linewidth=5, label = "SNpred")
# 	_S = LinRange(-0.25, 0.25, 1000)
# 	plot!([Hom.α(s)[1] for s in _S], [Hom.α(s)[2] for s in _S], linewidth=5, label = "HomPred")
#
# 	plot!(hopf_codim2, branchlabel = ["Hopf"], vars = (:β1, :β2), color = :black)
# 	plot!(br_hom, vars = (:β1, :param), branchlabel = ["Hom"], color=:yellow)
# 	xlims!(-0.02, 0.00001);ylims!(-0.1,0)
########################
br_hom = continuation(prob, btpt,
	# ShootingProblem(ODEProblem(Fbt!, zeros(2), (0., 1.), par), Rodas5P(), [zeros(2)]; abstol = 1e-12, reltol = 1e-11),
	ShootingProblem(ODEProblem(Fbt!, zeros(2), (0., 1.), par), Rodas5P(), [zeros(2) for _ in 1:5]; abstol = 1e-13, reltol = 1e-11, parallel = true),
	# PALC(tangent = Bordered()),
	PALC(),
	# MoorePenrose(),
	setproperties(optc_hom, maxSteps = 100, saveSolEveryStep = 1, dsmax = 9e-1, plotEveryStep = 1, detectEvent=0);
	amplitude = 5e-3, ϵ0 = 2e-4,
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	freeparams = ((@lens _.T), (@lens _.ϵ1)),
	# freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)), # WORK BEST
	# freeparams = ((@lens _.T),),
	verbosity = 0, plot = false,
	updateEveryStep = 1, # ne marche pas sinon. pb n est pas avec updatesection
	plotSolution = plotHom,
	callbackN = BK.cbMaxNorm(1e1),
	normC = norminf
	)

_sol = getHomoclinicOrbit(br_hom.prob.VF.F, br_hom.sol[end].x, BK.setParam(br_hom,  br_hom.sol[end].p))
	# plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[:,:]'))
