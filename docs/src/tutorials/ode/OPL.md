# Nonlinear laser model

```@contents
Pages = ["OPL.md"]
Depth = 3
```

The following model is taken from [^Pusuluri]:

$$\left\{\begin{aligned}
\dot{\beta} & =-\sigma \beta+g p_{23}, \\
\dot{p}_{21} & =-p_{21}-\beta p_{31}+a D_{21}, \\
\dot{p}_{23} & =-p_{23}+\beta D_{23}-a p_{31}, \\
\dot{p}_{31} & =-p_{31}+\beta p_{21}+a p_{23}, \\
\dot{D}_{21} & =-b\left(D_{21}-D_{21}^0\right)-4 a p_{21}-2 \beta p_{23}, \\
\dot{D}_{23} & =-b\left(D_{23}-D_{23}^0\right)-2 a p_{21}-4 \beta p_{23},
\end{aligned}\right.$$

It is easy to encode the ODE as follows

```@example TUTOPL
using Revise, Plots
using Setfield
using BifurcationKit
using HclinicBifurcationKit
const BK = BifurcationKit

record_from_solution(x, p) = (D₂₃ = x[6], β = x[1],)
####################################################################################################
function OPL!(dz, u, p, t = 0)
	(;b, σ, g, a, D₂₁⁰, D₂₃⁰)  = p
	β, p₂₁, p₂₃, p₃₁, D₂₁, D₂₃ = u
	dz[1] = -σ * β + g * p₂₃
	dz[2] =	-p₂₁ - β * p₃₁ + a * D₂₁
	dz[3] = -p₂₃ + β * D₂₃ - a * p₃₁
	dz[4] = -p₃₁ + β * p₂₁ + a * p₂₃
	dz[5] = -b * (D₂₁ - D₂₁⁰) - 4a * p₂₁ - 2β * p₂₃
	dz[6] = -b * (D₂₃ - D₂₃⁰) - 2a * p₂₁ - 4β * p₂₃
	dz
end

par_OPL = (b = 1.2, σ = 2.0, g=50., a = 1., D₂₁⁰ = -1., D₂₃⁰ = 0.)
z0 = zeros(6)
prob = BK.BifurcationProblem(OPL!, z0, par_OPL, (@lens _.a); record_from_solution)

nothing #hide
```

We first compute the branch of equilibria

```@example TUTOPL
opts_br = ContinuationPar(p_min = -1., p_max = 8., ds = 0.001, dsmax = 0.06, n_inversion = 6, nev = 6)
br = continuation(prob, PALC(), opts_br; normC = norminf)

plot(br, plotfold=true)

br2 = continuation(re_make(prob; u0 = [1.6931472491037485, -0.17634826359471437, 0.06772588996414994, -0.23085768742546342, -0.5672243219935907, -0.09634826359471438]), PALC(), opts_br; normC = norminf)
scene = plot(br, br2)
```

## Codimension 2 bifurcations


```@example TUTOPL
sn_br = continuation(br, 1, (@lens _.b), ContinuationPar(opts_br, detect_bifurcation = 1, max_steps = 80) ;
	detect_codim2_bifurcation = 2,
	start_with_eigen = true,
	bothside = true,
	)
show(sn_br)
```

```@example TUTOPL
hopf_br = continuation(br, 2, (@lens _.b), ContinuationPar(opts_br, detect_bifurcation = 1, max_steps = 140),
	detect_codim2_bifurcation = 2,
	bothside = true,
	)
show(hopf_br)
```

```@example TUTOPL
hopf_br2 = continuation(br2, 1, (@lens _.b), ContinuationPar(opts_br, detect_bifurcation = 1, max_steps = 140),
	detect_codim2_bifurcation = 2,
	bothside = true,
	)
show(hopf_br2)
```

```@example TUTOPL
plot(sn_br, vars = (:a, :b), branchlabel = "SN", )
plot!(hopf_br, branchlabel = "Hopf", vars = (:a, :b))
plot!(hopf_br2, branchlabel = "Hopf", vars = (:a, :b))
ylims!(0,1.5)
```

## Branch of homoclinic orbits with Orthogonal Collocation


```@example TUTOPL
# plotting function
function plotPO(x, p; k...)
	xtt = BK.get_periodic_orbit(p.prob, x, set(getparams(p.prob), BK.getlens(p.prob), p.p))
	plot!(xtt.t, xtt[1,:]; markersize = 2, k...)
	plot!(xtt.t, xtt[6,:]; k...)
	scatter!(xtt.t, xtt[1,:]; markersize = 1, legend = false, k...)
end

# record function
function recordPO(x, p)
	xtt = BK.get_periodic_orbit(p.prob, x, set(getparams(p.prob), BK.getlens(p.prob), p.p))
	period = BK.getperiod(p.prob, x, p.p)
	return (max = maximum(xtt[6,:]), min = minimum(xtt[6,:]), period = period, )
end

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  max_iterations = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, p_max = 6.8, p_min=-5., max_steps = 100, newton_options = (@set optn_po.tol = 1e-8), detect_bifurcation = 0, plot_every_step = 3)

br_coll = continuation(
	br2, 1,
	opts_po_cont,
	PeriodicOrbitOCollProblem(20, 4; meshadapt = true, update_section_every_step = 2);
	ampfactor = 1., δp = 0.0015,
	verbosity = 2,	plot = true,
	alg = PALC(tangent = Bordered()),
	record_from_solution = recordPO,
	plot_solution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		plot!(br,  subplot=1, putbifptlegend = false)
		plot!(br2, subplot=1, putbifptlegend = false)
		end,
	finalise_solution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
		# limit the period
			return z.u[end] < 150
			true
		end,
	normC = norminf)
	
_sol = get_periodic_orbit(br_coll, length(br_coll))
BK.plot(_sol.t, _sol.u'; marker = :d, markersize = 1,title = "Last periodic orbit on branch")
```

```@example TUTOPL
probhom, solh = generate_hom_problem(
    setproperties(br_coll.prob.prob, meshadapt=true, K = 100),
    br_coll.sol[end].x.sol,
    BK.setparam(br_coll, br_coll.sol[end].p),
    BK.getlens(br_coll);
    update_every_step = 4,
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

_sol = get_homoclinic_orbit(probhom, solh, BK.getparams(probhom);)
BK.plot(_sol.t, _sol.u'; marker = :d, markersize = 1,title = "Guess for homoclinic orbit")
```

```@example TUTOPL
optn_hom = NewtonPar(verbose = true, tol = 1e-10, max_iterations = 5)
optc_hom = ContinuationPar(newton_options = optn_hom, ds = -1e-4, dsmin = 1e-5, detect_bifurcation = 0, detect_event = 2, p_min = -1.01, dsmax = 4e-2, p_max = 1.5, max_steps = 100)

br_hom_c = continuation(
			deepcopy(probhom), solh, (@lens _.b),
			PALC(tangent = Bordered()),
			optc_hom;
	verbosity = 1, plot = true,
	normC = norminf,
	plot_solution = (x,p;k...) -> begin
		𝐇𝐨𝐦 = p.prob
		par0 = set(BK.getparams(𝐇𝐨𝐦), BK.getlens(𝐇𝐨𝐦), x.x[end][1])
		par0 = set(par0, (@lens _.b), p.p)
		sol = get_homoclinic_orbit(𝐇𝐨𝐦, x, par0)
		m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
		plot!(sol.t, sol[:,:]',subplot=3, markersize = 1, marker=m)
	end,
	)

plot(sn_br, vars = (:a, :b), branchlabel = "SN", )
plot!(hopf_br,  branchlabel = "AH₀",  vars = (:a, :b))
plot!(hopf_br2, branchlabel = "AH₁₂", vars = (:a, :b))
plot!(br_hom_c, branchlabel = "H₀",   vars = (:a, :b))
ylims!(0.1,1.5)
```


## Branch of homoclinic orbits with Multiple Shooting

```@example TUTOPL
using DifferentialEquations
probsh = ODEProblem(OPL!, copy(z0), (0., 1000.), par_OPL; abstol = 1e-12, reltol = 1e-10)

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  max_iterations = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.075, ds= -0.001, dsmin = 1e-4, p_max = 6.8, p_min=-5., max_steps = 130, newton_options = (@set optn_po.tol = 1e-8), tol_stability = 1e-4, detect_bifurcation = 0)

br_sh = continuation(
	br2, 1,
	opts_po_cont,
	ShootingProblem(8, probsh, Rodas5P(); parallel = true, abstol = 1e-13, reltol = 1e-11);
	ampfactor = 1., δp = 0.0015,
	verbosity = 2,	plot = true,
	record_from_solution = recordPO,
	callback_newton = BK.cbMaxNorm(1e0),
	plot_solution = (x, p; k...) -> begin
		plotPO(x, p; k...)
		## add previous branch
		# plot!(br, subplot=1, putbifptlegend = false)
		end,
	finalise_solution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
		# limit the period
			return z.u[end] < 100
			true
		end,
	normC = norminf)

_sol = get_periodic_orbit(br_sh.prob.prob, br_sh.sol[end].x, BK.setparam(br_sh,  br_sh.sol[end].p))
plot(_sol.t, _sol[1:3,:]')
```

```@example TUTOPL
# homoclinic
probhom, solh = generate_hom_problem(
	br_sh.prob.prob, br_sh.sol[end].x,
	BK.setparam(br_sh, br_sh.sol[end].p),
	BK.getlens(br_sh);
	verbose = true,
	update_every_step = 4,
	# ϵ0 = 1e-6, ϵ1 = 1e-5,
	t0 = 75, t1 = 25,
	# freeparams = ((@lens _.T), (@lens _.ϵ1),)
	# freeparams = ((@lens _.T), (@lens _.ϵ0)),
	# freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)), # WORK BEST
	freeparams = ((@lens _.T),),
	)

_sol = get_homoclinic_orbit(probhom, solh, BK.getparams(probhom); saveat=.1)
plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[1:4,:]'))

optn_hom = NewtonPar(verbose = true, tol = 1e-9, max_iterations = 7)
optc_hom = ContinuationPar(newton_options = optn_hom, ds = -1e-4, dsmin = 1e-6, max_steps = 300, detect_bifurcation = 0, dsmax = 12e-2, plot_every_step = 3, p_max = 7., detect_event = 2, a = 0.9)

br_hom_sh = continuation(
			deepcopy(probhom), solh, (@lens _.b),
			PALC(),
			optc_hom;
	verbosity = 3, plot = true,
	callback_newton = BK.cbMaxNorm(1e0),
	normC = norminf,
	plot_solution = (x,p;k...) -> begin
		𝐇𝐨𝐦 = p.prob
		par0 = set(BK.getparams(𝐇𝐨𝐦), BK.getlens(𝐇𝐨𝐦), x.x[end][1])
		par0 = set(par0, (@lens _.b), p.p)
		sol = get_homoclinic_orbit(𝐇𝐨𝐦, x, par0)
		m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
		plot!(sol.t, sol[1:6,:]',subplot=3, markersize = 1, marker=m)
	end,
	)

_sol = get_homoclinic_orbit(probhom, br_hom_sh.sol[end].x, BK.setparam(br_hom_sh, br_hom_sh.sol[end].p); saveat=.1)
plot(_sol.t, _sol[1:6,:]')
```

```@example TUTOPL
plot(sn_br, vars = (:a, :b), branchlabel = "SN", )
plot!(hopf_br, branchlabel = "AH₀", vars = (:a, :b))
plot!(hopf_br2, branchlabel = "AH₁₂", vars = (:a, :b))
plot!(br_hom_c, branchlabel = "Hc₀", vars = (:a, :b))
plot!(br_hom_sh, branchlabel = "Hsh₀", vars = (:a, :b), linewidth = 3)
ylims!(0,1.5)
```


## References

[^Pusuluri]:> Pusuluri, K, H G E Meijer, and A L Shilnikov. “Homoclinic Puzzles and Chaos in a Nonlinear Laser Model,” n.d.

