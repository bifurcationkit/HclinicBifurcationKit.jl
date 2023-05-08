# using Revise, Plots, AbbreviatedStackTraces
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

# plot(br, plotfold=true)
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
optn_po = NewtonPar(verbose = false, tol = 1e-8,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= -0.001, dsmin = 1e-4, pMax = 1.8, pMin=-5., maxSteps = 130, newtonOptions = (@set optn_po.tol = 1e-8), detectBifurcation = 0, plotEveryStep = 3, saveSolEveryStep=1,)

br_coll = continuation(
	br, 4, opts_po_cont,
	PeriodicOrbitOCollProblem(50, 4; meshadapt = true, updateSectionEveryStep = 2);
	ampfactor = 1., δp = 0.001,
	# verbosity = 2,	plot = true,
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

show(probhom)
HclinicBifurcationKit.generateHomoclinicSolution(probhom.bvp, t->ones(3), 1.)	
####################################################################################################
# same with shooting
using OrdinaryDiffEq
probsh = ODEProblem(freire!, copy(z0), (0., 1000.), par_freire; abstol = 1e-12, reltol = 1e-10)

# newton parameters
optn_po = NewtonPar(verbose = false, tol = 1e-8,  maxIter = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.15, ds= -0.0001, dsmin = 1e-4, pMax = 1.8, pMin=-5., maxSteps = 130, newtonOptions = (@set optn_po.tol = 1e-8), tolStability = 1e-4, detectBifurcation = 0, plotEveryStep = 20, saveSolEveryStep=1)

br_sh = continuation(
	br, 4, opts_po_cont,
	ShootingProblem(10, probsh, Rodas5P(); parallel = false);
	ampfactor = 1.0, δp = 0.001,
	# verbosity = 2,	plot = true,
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

# _sol = getPeriodicOrbit(br_sh.prob.prob, br_sh.sol[end].x, BK.setParam(br_sh,  br_sh.sol[end].p))
		# plot(_sol.t, _sol[:,:]')
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

show(probhom)