# Autonomous electronic circuit (aBS from BT)

```@contents
Pages = ["tutorialsFreire.md"]
Depth = 3
```

The following model is taken from [^Freire]:

$$\left\{\begin{aligned}
& r \dot{x}=-\nu x+\beta(y-x)-A_3 x^3+B_3(y-x)^3+\epsilon \\
& \dot{y}=-\beta(y-x)-z-B_3(y-x)^3 \\
& \dot{z}=y
\end{aligned}\right.$$

This is an example of computing the curve of homoclinics from a Bogdanov-Takens bifurcation point. 
It is easy to encode the ODE as follows

```@example TUTFREIRE
using Revise, Plots
using Parameters, Setfield, LinearAlgebra, Test, ForwardDiff
using BifurcationKit, Test
using HclinicBifurcationKit
const BK = BifurcationKit

norminf(x) = norm(x, Inf)
recordFromSolution(x, p) = (x = x[1], y = x[2])

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

nothing #hide
```

We first compute the branch of equilibria

```@example TUTFREIRE
opts_br = ContinuationPar(pMin = -1.4, pMax = 2.8, ds = 0.001, dsmax = 0.05, nInversion = 6, detectBifurcation = 3, maxBisectionSteps = 25, nev = 3)

br = continuation(prob, PALC(tangent = Bordered()), opts_br, normC = norminf)

scene = plot(br, plotfold=true)
```

With detailed information:

```@example TUTFREIRE
br
```

## Codimension 2 bifurcations


```@example TUTFREIRE
sn_br = continuation(br, 2, (@lens _.ν), ContinuationPar(opts_br, detectBifurcation = 1, saveSolEveryStep = 1, dsmax = 0.01, maxSteps = 80) ;
	alg = PALC(),
	detectCodim2Bifurcation = 2,
	startWithEigen = true,
	updateMinAugEveryStep = 1,
	bothside = true,
	)

hopf_br = continuation(br, 4, (@lens _.ν), ContinuationPar(opts_br, detectBifurcation = 1, saveSolEveryStep = 1, maxSteps = 140, dsmax = 0.02, nInversion = 6),
	detectCodim2Bifurcation = 2,
	startWithEigen = true,
	updateMinAugEveryStep = 1,
	bothside = true,
	)

plot(sn_br, hopf_br, ylims = (-0.1, 1.25))
```

## Branch of homoclinic orbits with Orthogonal Collocation

Plotting function:

```@example TUTFREIRE
function plotHom(x,p;k...)
	𝐇𝐨𝐦 = p.prob
	par0 = set(BK.getParams(𝐇𝐨𝐦), BK.getLens(𝐇𝐨𝐦), x.x[end][1])
	par0 = set(par0, p.lens, p.p)
	sol = getHomoclinicOrbit(𝐇𝐨𝐦, x, par0)
	m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
	plot!(sol.t, sol[:,:]',subplot=3, markersize = 1, marker=m)
end
```

Branching from BT point

```@example TUTFREIRE
# Bogdanov-Takens bifurcation point
btpt = getNormalForm(sn_br, 2; nev = 3, autodiff = false)

# curve of homoclinics
br_hom_c = continuation(
			prob,
			btpt,
			# we use mesh adaptation
			PeriodicOrbitOCollProblem(50, 3; meshadapt = true, K = 100),
			PALC(tangent = Bordered()),
			setproperties(opts_br, maxSteps = 8, saveSolEveryStep = 1, dsmax = 1e-2, plotEveryStep = 1, pMin = -1.01, ds = 0.01, detectEvent = 2, detectBifurcation = 0);
	verbosity = 1, plot = true,
	ϵ0 = 1e-5, amplitude = 2e-3,
	# freeparams = ((@lens _.T), (@lens _.ϵ1),)
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)),
	normC = norminf,
	plotSolution = plotHom,
	updateEveryStep = 4,
	)
title!("")	
```

```@example TUTFREIRE
br_hom_c = continuation(
			prob,
			btpt,
			# we use mesh adaptation
			PeriodicOrbitOCollProblem(50, 3; meshadapt = true, K = 100),
			PALC(tangent = Bordered()),
			setproperties(opts_br, maxSteps = 130, saveSolEveryStep = 1, dsmax = 1e-2, plotEveryStep = 1, pMin = -1.01, ds = 0.01, detectEvent = 2, detectBifurcation = 0);
	verbosity = 0, plot = false,
	ϵ0 = 1e-5, amplitude = 2e-3,
	# freeparams = ((@lens _.T), (@lens _.ϵ1),)
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)),
	normC = norminf,
	updateEveryStep = 4,
	)
plot(sn_br, hopf_br, ylims = (0, 1.25), branchlabel = ["SN", "HOPF"])
	plot!(br_hom_c, branchlabel = "Hom")
```

## Branch of homoclinic orbits with Multiple Shooting

```@example TUTFREIRE
using DifferentialEquations
probsh = ODEProblem(freire!, copy(z0), (0., 1000.), par_freire; abstol = 1e-12, reltol = 1e-10)

optn_hom = NewtonPar(verbose = true, tol = 1e-10, maxIter = 7)
optc_hom = ContinuationPar(newtonOptions = optn_hom, ds = 1e-4, dsmin = 1e-6, dsmax = 1e-3, plotEveryStep = 1,maxSteps = 10, detectBifurcation = 0, saveSolEveryStep = 1)

br_hom_sh = continuation(
			prob,
			btpt,
			ShootingProblem(12, probsh, Rodas5P(); parallel = true, abstol = 1e-13, reltol = 1e-12),
			PALC(tangent = Bordered()),
			setproperties(optc_hom, maxSteps = 200, saveSolEveryStep = 1, ds = 1e-3, dsmax = 3e-2, plotEveryStep = 50, detectEvent = 2, a = 0.9, pMin = -1.01);
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
title!("")
```


## References

[^Freire] :> Freire, E., A.J. Rodríguez-Luis, E. Gamero, and E. Ponce. “A Case Study for Homoclinic Chaos in an Autonomous Electronic Circuit.” Physica D: Nonlinear Phenomena 62, no. 1–4 (January 1993): 230–53. https://doi.org/10.1016/0167-2789(93)90284-8.
