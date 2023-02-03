using Revise, Plots
using Parameters, Setfield, LinearAlgebra, Test, ForwardDiff
using BifurcationKit, Test
using HclinicBifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)
recordFromSolution(x, p) = (x = x[1], y = x[2])
####################################################################################################
function freire!(dz, u, p, t)
	@unpack ν, β, A₃, B₃, r, ϵ = p
	x, y, z = u
	dz[1] = (-ν*x + β*(y-x) - A₃*x^3 + B₃*(y-x)^3 + ϵ)/r
	dz[2] =	-β*(y-x) - z - B₃*(y-x)^3
	dz[3] = y
	dz
end

freire(z, p) = freire!(similar(z), z, p, 0)
par_freire = (ν = -0.75, β = -0.1, A₃ = 0.328578, B₃ = 0.933578, r = 0.6, ϵ = 0.01)
z0 = [0.7,0.3,0.1]
z0 = zeros(3)
prob = BK.BifurcationProblem(freire, z0, par_freire, (@lens _.β); recordFromSolution = recordFromSolution)

opts_br = ContinuationPar(pMin = -1.4, pMax = 2.8, ds = 0.001, dsmax = 0.05, nInversion = 6, detectBifurcation = 3, maxBisectionSteps = 25, nev = 3, maxSteps = 2000)
	@set! opts_br.newtonOptions.verbose = false
	br = continuation(prob, PALC(tangent = Bordered()), opts_br; verbosity = 0,
	bothside = false, normC = norminf)

plot(br, plotfold=true)
####################################################################################################
# DefCon
alg = BK.DefCont(deflationOperator = DeflationOperator(2, dot, .001, [z0]), perturbSolution = (x,p,id) -> (x  .+ 0.03 .* rand(length(x))))

opts_br_dc = ContinuationPar(pMin = -1.4, pMax = 2.8, ds = 0.001, dsmax = 0.01, detectBifurcation = 3, maxBisectionSteps = 25, nev = 3, plotEveryStep = 100, maxSteps = 1000)
	@set! opts_br_dc.newtonOptions.verbose = false
	brdc = continuation(prob, alg, opts_br_dc; verbosity = 0,
	normN = norminf)
####################################################################################################
sn_br = continuation(br, 2, (@lens _.ν), ContinuationPar(opts_br, detectBifurcation = 1, saveSolEveryStep = 1, dsmax = 0.01, maxSteps = 80) ;
	alg = PALC(),
	detectCodim2Bifurcation = 2,
	startWithEigen = true,
	updateMinAugEveryStep = 1,
	bothside = true,
	)

plot(sn_br)
bt = getNormalForm(sn_br, 2, verbose = true, detailed = true, autodiff = false)

hopf_br = continuation(br, 4, (@lens _.ν), ContinuationPar(opts_br, detectBifurcation = 1, saveSolEveryStep = 1, maxSteps = 140, dsmax = 0.02, nInversion = 6),
	detectCodim2Bifurcation = 2,
	startWithEigen = true,
	updateMinAugEveryStep = 1,
	bothside = true,
	)

plot(sn_br, hopf_br, ylims = (-0.1, 1.25))
####################################################################################################
function plotHom(x,p;k...)
	𝐇𝐨𝐦 = p.prob
	par0 = set(BK.getParams(𝐇𝐨𝐦), BK.getLens(𝐇𝐨𝐦), x.x[end][1])
	par0 = set(par0, p.lens, p.p)
	sol = getHomoclinicOrbit(𝐇𝐨𝐦, x, par0)
	m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
	plot!(sol.t, sol[:,:]',subplot=3, markersize = 1, marker=m)
end

btpt = getNormalForm(sn_br, 2; nev = 3, autodiff = false)

br_hom_c = continuation(
			prob,
			btpt,
			PeriodicOrbitOCollProblem(50, 3; meshadapt = true, K = 100),
			PALC(tangent = Bordered()),
			setproperties(opts_br, maxSteps = 130, saveSolEveryStep = 1, dsmax = 1e-2, plotEveryStep = 1, pMin = -1.01, ds = 0.01, detectEvent = 0, detectBifurcation = 0);
	verbosity = 1, plot = true,
	ϵ0 = 1e-5, amplitude = 2e-3,
	# freeparams = ((@lens _.T), (@lens _.ϵ1),)
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)),
	normC = norminf,
	plotSolution = plotHom,
	updateEveryStep = 4,
	)

plot(sn_br, hopf_br, ylims = (0, 1.25), branchlabel = ["SN", "HOPF"])
	plot!(br_hom_c, branchlabel = "Hom")
####################################################################################################
# plotting function
function plotPO(x, p; k...)
	xtt = BK.getPeriodicOrbit(p.prob, x, @set par_freire.β = p.p)
	plot!(xtt.t, xtt[1,:]; markersize = 2, k...)
	plot!(xtt.t, xtt[2,:]; k...)
	plot!(xtt.t, xtt[3,:]; marker = :d, markersize = 1, legend = false, k...)
end

# record function
function recordPO(x, p)
	xtt = BK.getPeriodicOrbit(p.prob, x, @set par_freire.β = p.p)
	period = BK.getPeriod(p.prob, x, p.p)
	return (period = period, max = maximum(xtt[1,:]), min = minimum(xtt[1,:]))
end

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= -0.001, dsmin = 1e-4, pMax = 1.8, pMin=-5., maxSteps = 130, newtonOptions = (@set optn_po.tol = 1e-8), detectBifurcation = 0, plotEveryStep = 3, saveSolEveryStep=1,)

br_coll = continuation(
	br, 4, opts_po_cont,
	PeriodicOrbitOCollProblem(50, 4; meshadapt = true, updateSectionEveryStep = 2);
	ampfactor = 1., δp = 0.01,
	verbosity = 2,	plot = true,
	alg = PALC(tangent = Bordered()),
	recordFromSolution = recordPO,
	plotSolution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		# plot!(br, subplot=1, putbifptlegend = false)
		end,
	finaliseSolution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
		# limit the period
			return z.u[end] < 200
			true
		end,
	normC = norminf)

_sol = getPeriodicOrbit(br_coll.prob.prob, br_coll.sol[end].x,0)
		plot(_sol.t, _sol[:,:]')
####################################################################################################
# homoclinic
probhom, solh = generateHomProblem(
	setproperties(br_coll.prob.prob, meshadapt=true, K = 100),
	br_coll.sol[end].x,
	BK.setParam(br_coll, br_coll.sol[end].p),
	BK.getLens(br_coll);
	updateEveryStep = 4,
	ϵ0 = 1e-9,
	ϵ1 = 1e-6,
	# freeparams = ((@lens _.T), (@lens _.ϵ1),)
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)),
	# freeparams = ((@lens _.T),),
	)

#####

_sol = getHomoclinicOrbit(probhom, solh, BK.getParams(probhom);)
	plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[:,:]'))

optn_hom = NewtonPar(verbose = true, tol = 1e-10, maxIter = 5)
optc_hom = ContinuationPar(newtonOptions = optn_hom, ds = 0.0001, dsmin = 1e-5, plotEveryStep = 10, maxSteps = 100, detectBifurcation = 0, detectEvent = 2, saveSolEveryStep = 1, pMin = -1.01)

solh.x[2] .=0

br_hom_c = continuation(
			deepcopy(probhom), solh, (@lens _.ν),
			PALC(tangent = Bordered()),
			setproperties(optc_hom, maxSteps = 130, saveSolEveryStep = 1, dsmax = 1e-2, plotEveryStep = 1);
	verbosity = 4, plot = true,
	normC = norminf,
	plotSolution = plotHom,
	)

using PrettyTables
br_hom_c.branch[end-20:end] |> pretty_table

plot(sn_br, hopf_br, ylims = (0, 1.25))
	plot!(br_hom_c)

_sol = getHomoclinicOrbit(probhom, br_hom_c.sol[end].x, BK.setParam(br_hom_c,  br_hom_c.sol[end].p))
	plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[:,:]'))
####################################################################################################
# same with shooting
using DifferentialEquations
probsh = ODEProblem(freire!, copy(z0), (0., 1000.), par_freire; abstol = 1e-12, reltol = 1e-10)

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.15, ds= -0.0001, dsmin = 1e-4, pMax = 1.8, pMin=-5., maxSteps = 130, newtonOptions = (@set optn_po.tol = 1e-8), tolStability = 1e-4, detectBifurcation = 0, plotEveryStep = 20, saveSolEveryStep=1)

br_sh = continuation(
	br, 4, opts_po_cont,
	ShootingProblem(10, probsh, Rodas5P(); parallel = true, abstol = 1e-13, reltol = 1e-12);
	ampfactor = 1., δp = 0.01,
	verbosity = 2,	plot = true,
	recordFromSolution = recordPO,
	plotSolution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		# plot!(br, subplot=1, putbifptlegend = false)
		end,
	finaliseSolution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
		# limit the period
			return z.u[end] < 300
			true
		end,
	normC = norminf)

_sol = getPeriodicOrbit(br_sh.prob.prob, br_sh.sol[end].x, BK.setParam(br_sh,  br_sh.sol[end].p))
		plot(_sol.t, _sol[:,:]')
#######################################
# homoclinic
probhom, solh = generateHomProblem(
	br_sh.prob.prob, br_sh.sol[end].x,
	BK.setParam(br_sh, br_sh.sol[end].p),
	BK.getLens(br_sh);
	verbose = true,
	updateEveryStep = 4,
	ϵ0 = 7e-8,
	ϵ1 = 8e-8,
	# freeparams = ((@lens _.T), (@lens _.ϵ1),)
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)), # WORK BEST
	# freeparams = ((@lens _.T),),
	)

solh.x[2] .=0

_sol = getHomoclinicOrbit(probhom, solh, BK.getParams(probhom); saveat=.1)
	plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[:,:]'))

optn_hom = NewtonPar(verbose = true, tol = 1e-10, maxIter = 7)
optc_hom = ContinuationPar(newtonOptions = optn_hom, ds = 1e-4, dsmin = 1e-6, dsmax = 1e-3, plotEveryStep = 1,maxSteps = 10, detectBifurcation = 0, saveSolEveryStep = 1)

br_hom_sh = continuation(
			deepcopy(probhom), solh, (@lens _.ν),
			PALC(tangent = Bordered()),
			# MoorePenrose(),
			setproperties(optc_hom, maxSteps = 20380, saveSolEveryStep = 1, dsmax = 3e-2, plotEveryStep = 10, detectEvent = 2);
	verbosity = 4, plot = true,
	callbackN = BK.cbMaxNorm(1e1),
	normC = norminf,
	plotSolution = plotHom,
	)

_sol = getHomoclinicOrbit(br_hom_sh.prob.VF.F, br_hom_sh.sol[end].x, BK.setParam(br_hom_sh, br_hom_sh.sol[end].p); saveat=.1)
	plot(_sol.t, _sol[:,:]')

plot(hopf_br, br_hom_c, br_hom_sh)

#######
br_hom_sh = continuation(
			prob,
			btpt,
			ShootingProblem(12, probsh, Rodas5P(); parallel = true, abstol = 1e-13, reltol = 1e-12),
			PALC(tangent = Bordered()),
			setproperties(optc_hom, maxSteps = 15000, saveSolEveryStep = 1, ds = 1e-3, dsmax = 3e-2, plotEveryStep = 50, detectEvent = 2, a = 0.9, pMin = -1.01);
	verbosity = 1, plot = true,
	ϵ0 = 1e-6, amplitude = 2e-2,
	updateEveryStep = 2,
	# freeparams = ((@lens _.T), (@lens _.ϵ1),),
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	# freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)),
	freeparams = ((@lens _.T),),
	normC = norminf,
	plotSolution = plotHom,
	maxT = 45.,
	callbackN = BK.cbMaxNorm(1e1),
	)
