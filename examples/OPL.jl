using Revise, Plots
using Parameters, Setfield, LinearAlgebra, Test, ForwardDiff
using BifurcationKit, Test
using HclinicBifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)
recordFromSolution(x, p) = (D₂₃ = x[6], β = x[1],)
####################################################################################################
function OPL!(dz, u, p, t)
	@unpack b, σ, g, a, D₂₁⁰, D₂₃⁰  = p
	β, p₂₁, p₂₃, p₃₁, D₂₁, D₂₃ = u
	dz[1] = -σ * β + g * p₂₃
	dz[2] =	-p₂₁ - β * p₃₁ + a * D₂₁
	dz[3] = -p₂₃ + β * D₂₃ - a * p₃₁
	dz[4] = -p₃₁ + β * p₂₁ + a * p₂₃
	dz[5] = -b * (D₂₁ - D₂₁⁰) - 4a * p₂₁ - 2β * p₂₃
	dz[6] = -b * (D₂₃ - D₂₃⁰) - 2a * p₂₁ - 4β * p₂₃
	dz
end

OPL(z, p) = OPL!(similar(z), z, p, 0)
par_OPL = (b = 1.2, σ = 2.0, g=50., a = 1., D₂₁⁰ = -1., D₂₃⁰ = 0.)
z0 = zeros(6)
prob = BK.BifurcationProblem(OPL, z0, par_OPL, (@lens _.a); recordFromSolution = recordFromSolution)

opts_br = ContinuationPar(pMin = -1., pMax = 8., ds = 0.001, dsmax = 0.06, nInversion = 6, detectBifurcation = 3, maxBisectionSteps = 25, nev = 6, plotEveryStep = 20, maxSteps = 100, saveSolEveryStep = 1, detectFold = true)
	opts_br = @set opts_br.newtonOptions.verbose = false
	br = continuation(prob, PALC(tangent = Secant()), opts_br;
	bothside = false, normC = norminf)

plot(br, plotfold=true)

br2 = continuation(reMake(prob; u0 = [1.6931472491037485, -0.17634826359471437, 0.06772588996414994, -0.23085768742546342, -0.5672243219935907, -0.09634826359471438]), PALC(tangent = Secant()), opts_br; verbosity = 0, bothside = false, normC = norminf)
plot(br, br2)
####################################################################################################
sn_br = continuation(br, 1, (@lens _.b), ContinuationPar(opts_br, detectBifurcation = 1, saveSolEveryStep = 1, maxSteps = 80) ;
	alg = PALC(),
	verbosity = 0,
	detectCodim2Bifurcation = 2,
	startWithEigen = true,
	updateMinAugEveryStep = 1,
	bothside = true,
	)

plot(sn_br)
bt = getNormalForm(sn_br, 2, verbose = true, detailed = true, autodiff = false)

hopf_br = continuation(br, 2, (@lens _.b), ContinuationPar(opts_br, detectBifurcation = 1, saveSolEveryStep = 1, maxSteps = 140),
	detectCodim2Bifurcation = 2,
	startWithEigen = true,
	updateMinAugEveryStep = 1,
	bothside = true,
	)

hopf_br2 = continuation(br2, 1, (@lens _.b), ContinuationPar(opts_br, detectBifurcation = 1, saveSolEveryStep = 1, maxSteps = 140),
	detectCodim2Bifurcation = 2,
	startWithEigen = true,
	updateMinAugEveryStep = 1,
	bothside = true,
	)

plot(sn_br, vars = (:a, :b), branchlabel = ["SN"], )
	plot!(hopf_br, branchlabel = ["Hopf"], vars = (:a, :b))
	plot!(hopf_br2, branchlabel = ["Hopf"], vars = (:a, :b))
	ylims!(0,1.5)

####################################################################################################
# plotting function
function plotPO(x, p; k...)
	xtt = BK.getPeriodicOrbit(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p))
	plot!(xtt.t, xtt[1,:]; markersize = 2, k...)
	plot!(xtt.t, xtt[6,:]; k...)
	scatter!(xtt.t, xtt[1,:]; markersize = 1, legend = false, k...)
end

# record function
function recordPO(x, p)
	xtt = BK.getPeriodicOrbit(p.prob, x, set(getParams(p.prob), BK.getLens(p.prob), p.p))
	period = BK.getPeriod(p.prob, x, p.p)
	return (max = maximum(xtt[6,:]), min = minimum(xtt[6,:]), period = period, )
end

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, pMax = 6.8, pMin=-5., maxSteps = 100, newtonOptions = (@set optn_po.tol = 1e-8), detectBifurcation = 0, plotEveryStep = 3, saveSolEveryStep=1,)

br_coll = continuation(
	# br, 2,
	br2, 1,
	opts_po_cont,
	PeriodicOrbitOCollProblem(20, 4; meshadapt = true, updateSectionEveryStep = 2);
	ampfactor = 1., δp = 0.0015,
	verbosity = 2,	plot = true,
	alg = PALC(tangent = Bordered()),
	recordFromSolution = recordPO,
	plotSolution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		plot!(br, subplot=1, putbifptlegend = false)
		plot!(br2, subplot=1, putbifptlegend = false)
		end,
	finaliseSolution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
		# limit the period
			return z.u[end] < 150
			true
		end,
	normC = norminf)

_sol = getPeriodicOrbit(br_coll.prob.prob, br_coll.sol[end].x, 0)
		plot(_sol.t, _sol[:,:]', marker = :d, markersize = 1)

####################################################################################################
# homoclinic
probhom, solh = generateHomProblem(
	setproperties(br_coll.prob.prob, meshadapt=true, K = 100),
	br_coll.sol[end].x,
	BK.setParam(br_coll, br_coll.sol[end].p),
	BK.getLens(br_coll);
	updateEveryStep = 4,
	verbose = true,
	# ϵ0 = 1e-7, ϵ1 = 1e-7, # WORK BEST
	# ϵ0 = 1e-8, ϵ1 = 1e-5, # maxT = 70,
	t0 = 0., t1 = 120.,
	# freeparams = ((@lens _.T), (@lens _.ϵ1),)
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)), # WORK BEST
	# freeparams = ((@lens _.T),),
	testOrbitFlip = false,
	testInclinationFlip = false
	)

#####

	_sol = getHomoclinicOrbit(probhom, solh, BK.getParams(probhom);)
	plot(_sol.t, _sol[:,:]', marker = :d, markersize = 1)

optn_hom = NewtonPar(verbose = true, tol = 1e-10, maxIter = 5)
optc_hom = ContinuationPar(newtonOptions = optn_hom, ds = -0.0001, dsmin = 1e-5, plotEveryStep = 10, maxSteps = 100, detectBifurcation = 0, detectEvent = 2, saveSolEveryStep = 1, pMin = -1.01)

br_hom_c = continuation(
			deepcopy(probhom), solh, (@lens _.b),
			PALC(tangent = Bordered()),
			# PALC(),
			# MoorePenrose(),
			setproperties(optc_hom, maxSteps = 100, saveSolEveryStep = 1, dsmax = 4e-2, plotEveryStep = 10, pMax = 1.5);
	verbosity = 1, plot = true,
	# callbackN = BK.cbMaxNorm(1e1),
	# bothside = true,
	normC = norminf,
	plotSolution = (x,p;k...) -> begin
		𝐇𝐨𝐦 = p.prob
		par0 = set(BK.getParams(𝐇𝐨𝐦), BK.getLens(𝐇𝐨𝐦), x.x[end][1])
		par0 = set(par0, (@lens _.b), p.p)
		sol = getHomoclinicOrbit(𝐇𝐨𝐦, x, par0)
		m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
		plot!(sol.t, sol[:,:]',subplot=3, markersize = 1, marker=m)
	end,
	)

using PrettyTables
br_hom_c.branch |> pretty_table

plot(sn_br, vars = (:a, :b), branchlabel = "SN", )
	plot!(hopf_br, branchlabel = "AH₀", vars = (:a, :b))
	plot!(hopf_br2, branchlabel = "AH₁₂", vars = (:a, :b))
	plot!(br_hom_c, branchlabel = "H₀", vars = (:a, :b))
	ylims!(0,1.5)
####################################################################################################
# same with shooting
using DifferentialEquations
probsh = ODEProblem(OPL!, copy(z0), (0., 1000.), par_OPL; abstol = 1e-12, reltol = 1e-10)

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.075, ds= -0.001, dsmin = 1e-4, pMax = 6.8, pMin=-5., maxSteps = 130, newtonOptions = (@set optn_po.tol = 1e-8), tolStability = 1e-4, detectBifurcation = 0, plotEveryStep = 10, saveSolEveryStep=1)

br_sh = continuation(
	# br, 2,
	br2, 1,
	opts_po_cont,
	ShootingProblem(8, probsh, Rodas5P(); parallel = true, abstol = 1e-13, reltol = 1e-11);
	ampfactor = 1., δp = 0.0015,
	verbosity = 2,	plot = true,
	recordFromSolution = recordPO,
	# alg = MoorePenrose(),
	callbackN = BK.cbMaxNorm(1e0),
	plotSolution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		# plot!(br, subplot=1, putbifptlegend = false)
		end,
	finaliseSolution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
		# limit the period
			return z.u[end] < 100
			true
		end,
	normC = norminf)

_sol = getPeriodicOrbit(br_sh.prob.prob, br_sh.sol[end].x, br_sh.sol[end].p)
		plot(_sol.t, _sol[:,:]')

#######################################
# homoclinic
probhom, solh = generateHomProblem(
	br_sh.prob.prob, br_sh.sol[end].x,
	BK.setParam(br_sh, br_sh.sol[end].p),
	BK.getLens(br_sh);
	verbose = true,
	updateEveryStep = 4,
	# ϵ0 = 1e-6, ϵ1 = 1e-5,
	t0 = 75, t1 = 25,
	# freeparams = ((@lens _.T), (@lens _.ϵ1),)
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	# freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)), # WORK BEST
	freeparams = ((@lens _.T),),
	)

	_sol = getHomoclinicOrbit(probhom, solh, BK.getParams(probhom); saveat=.1)
	plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[:,:]'))

optn_hom = NewtonPar(verbose = true, tol = 1e-9, maxIter = 7)
	optc_hom = ContinuationPar(newtonOptions = optn_hom, ds = -1e-4, dsmin = 1e-6, dsmax = 1e-3, plotEveryStep = 1, maxSteps = 10, detectBifurcation = 0, saveSolEveryStep = 1)

br_hom_sh = continuation(
			deepcopy(probhom), solh, (@lens _.b),
			# PALC(tangent = Bordered()),
			PALC(),
			# ANM(6, 1e-8)
			# MoorePenrose(),
			setproperties(optc_hom, maxSteps = 600, saveSolEveryStep = 1, dsmax = 12e-2, plotEveryStep = 3, pMax = 7., detectEvent = 2, a = 0.9);
	verbosity = 3, plot = true,
	callbackN = BK.cbMaxNorm(1e0),
	normC = norminf,
	plotSolution = (x,p;k...) -> begin
		𝐇𝐨𝐦 = p.prob
		par0 = set(BK.getParams(𝐇𝐨𝐦), BK.getLens(𝐇𝐨𝐦), x.x[end][1])
		par0 = set(par0, (@lens _.b), p.p)
		sol = getHomoclinicOrbit(𝐇𝐨𝐦, x, par0)
		m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
		plot!(sol.t, sol[:,:]',subplot=3, markersize = 1, marker=m)
	end,
	)

_sol = getHomoclinicOrbit(probhom, br_hom_sh.sol[end].x, BK.setParam(br_hom_sh, br_hom_sh.sol[end].p); saveat=.1)
	plot(_sol.t, _sol[:,:]')


plot(sn_br, vars = (:a, :b), branchlabel = "SN", )
	plot!(hopf_br, branchlabel = "AH₀", vars = (:a, :b))
	plot!(hopf_br2, branchlabel = "AH₁₂", vars = (:a, :b))
	plot!(br_hom_c, branchlabel = "Hc₀", vars = (:a, :b))
	plot!(br_hom_sh, branchlabel = "Hsh₀", vars = (:a, :b), linewidth = 3)
	ylims!(0,1.5)
