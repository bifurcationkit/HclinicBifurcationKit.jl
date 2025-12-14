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
using LinearAlgebra, Test, ForwardDiff
using BifurcationKit, Test
using HclinicBifurcationKit
const BK = BifurcationKit

recordFromSolution(x, p; k...) = (x = x[1], y = x[2])

function freire!(dz, u, p, t = 0)
	(;ν, β, A₃, B₃, r, ϵ) = p
	x, y, z = u
	dz[1] = (-ν*x + β*(y-x) - A₃*x^3 + B₃*(y-x)^3 + ϵ)/r
	dz[2] =	-β*(y-x) - z - B₃*(y-x)^3
	dz[3] = y
	dz
end

par_freire = (ν = -0.75, β = -0.1, A₃ = 0.328578, B₃ = 0.933578, r = 0.6, ϵ = 0.01)
z0 = [0.7,0.3,0.1]
z0 = zeros(3)
prob = BK.BifurcationProblem(freire!, z0, par_freire, (@optic _.β); record_from_solution = recordFromSolution)

nothing #hide
```

We first compute the branch of equilibria

```@example TUTFREIRE
opts_br = ContinuationPar(p_min = -1.4, p_max = 2.8, ds = 0.001, dsmax = 0.05, n_inversion = 6, nev = 3)

br = continuation(prob, PALC(tangent = Bordered()), opts_br, normC = norminf)

scene = plot(br, plotfold = true)
```

With detailed information:

```@example TUTFREIRE
br
```

## Codimension 2 bifurcations


```@example TUTFREIRE
sn_br = continuation(br, 2, (@optic _.ν), ContinuationPar(opts_br, detect_bifurcation = 1, dsmax = 0.01, max_steps = 80) ;
	alg = PALC(),
	detect_codim2_bifurcation = 2,
	start_with_eigen = true,
	update_minaug_every_step = 1,
	bothside = true,
	)

hopf_br = continuation(br, 4, (@optic _.ν), ContinuationPar(opts_br, detect_bifurcation = 1, save_sol_every_step = 1, max_steps = 140, dsmax = 0.02, n_inversion = 6),
	detect_codim2_bifurcation = 2,
	start_with_eigen = true,
	update_minaug_every_step = 1,
	bothside = true,
	)

plot(sn_br, hopf_br, ylims = (-0.1, 1.25))
```

## Branch of homoclinic orbits with Orthogonal Collocation

Plotting function:

```@example TUTFREIRE
function plotHom(x,p;k...)
	𝐇𝐨𝐦 = p.prob
	par0 = set(BK.getparams(𝐇𝐨𝐦), BK.getlens(𝐇𝐨𝐦), x.x[end][1])
	par0 = set(par0, p.lens, p.p)
	sol = get_homoclinic_orbit(𝐇𝐨𝐦, x, par0)
	m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
	plot!(sol.t, sol[1:3,:]',subplot=3, markersize = 1, marker=m)
end
```

Branching from BT point

```@example TUTFREIRE
# Bogdanov-Takens bifurcation point
btpt = get_normal_form(sn_br, 2; nev = 3, autodiff = false)

# curve of homoclinics
br_hom_c = continuation(
			prob,
			btpt,
			# we use mesh adaptation
			PeriodicOrbitOCollProblem(50, 3; meshadapt = false, K = 200),
			PALC(tangent = Bordered()),
			setproperties(opts_br, max_steps = 30, dsmax = 1e-2, plot_every_step = 1, p_min = -1.01, ds = 0.001, detect_event = 2, detect_bifurcation = 0);
	verbosity = 0, plot = false,
	ϵ0 = 1e-5, amplitude = 2e-3,
	# freeparams = ((@optic _.T), (@optic _.ϵ1),)
	# freeparams = ((@optic _.T), (@optic _.ϵ0)),
	freeparams = ((@optic _.ϵ0), (@optic _.ϵ1)),
	normC = norminf,
	update_every_step = 4,
	)
plot(sn_br, hopf_br, ylims = (0, 1.25), branchlabel = ["SN", "HOPF"])
plot!(br_hom_c, branchlabel = "Hom")
```

## Branch of homoclinic orbits with Multiple Shooting

```@example TUTFREIRE
using OrdinaryDiffEq
probsh = ODEProblem(freire!, copy(z0), (0., 1000.), par_freire; abstol = 1e-12, reltol = 1e-10)

optn_hom = NewtonPar(verbose = true, tol = 1e-10, max_iterations = 7)
optc_hom = ContinuationPar(newton_options = optn_hom, ds = 1e-3, dsmin = 1e-6, dsmax = 3e-2, plot_every_step = 50, max_steps = 200, detect_bifurcation = 0)

br_hom_sh = continuation(
			prob,
			btpt,
			ShootingProblem(12, probsh, Rodas5P(); parallel = true, abstol = 1e-13, reltol = 1e-12),
			PALC(tangent = Bordered()),
			setproperties(optc_hom, detect_event = 2, a = 0.9, p_min = -1.01);
	verbosity = 1, plot = true,
	ϵ0 = 1e-6, amplitude = 2e-2,
	update_every_step = 2,
	# freeparams = ((@optic _.T), (@optic _.ϵ1),),
	# freeparams = ((@optic _.T), (@optic _.ϵ0)),
	# freeparams = ((@optic _.ϵ0), (@optic _.ϵ1)),
	freeparams = ((@optic _.T),),
	normC = norminf,
	plot_solution = plotHom,
	maxT = 45.,
	callback_newton = BK.cbMaxNorm(1e1),
	)
title!("")
```


## References

[^Freire] :> Freire, E., A.J. Rodríguez-Luis, E. Gamero, and E. Ponce. “A Case Study for Homoclinic Chaos in an Autonomous Electronic Circuit.” Physica D: Nonlinear Phenomena 62, no. 1–4 (January 1993): 230–53. https://doi.org/10.1016/0167-2789(93)90284-8.
