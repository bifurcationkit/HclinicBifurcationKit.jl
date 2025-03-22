# using Revise, Plots, AbbreviatedStackTraces
using LinearAlgebra, Test, ForwardDiff
using BifurcationKit, Test
using HclinicBifurcationKit
const BK = BifurcationKit

recordFromSolution(x, p; k...) = (x = x[1], y = x[2])
####################################################################################################
function freire!(dz, u, p, t = 0)
    (;ν, β, A₃, B₃, r, ϵ) = p
    x, y, z = u
    dz[1] = (-ν*x + β*(y-x) - A₃*x^3 + B₃*(y-x)^3 + ϵ)/r
    dz[2] =    -β*(y-x) - z - B₃*(y-x)^3
    dz[3] = y
    dz
end

par_freire = (ν = -0.75, β = -0.1, A₃ = 0.328578, B₃ = 0.933578, r = 0.6, ϵ = 0.01)
z0 = [0.7,0.3,0.1]
z0 = zeros(3)
prob = BK.BifurcationProblem(freire!, z0, par_freire, (@optic _.β); record_from_solution = recordFromSolution)

opts_br = ContinuationPar(p_min = -1.4, p_max = 2.8, ds = 0.001, dsmax = 0.05, n_inversion = 6, detect_bifurcation = 3, max_bisection_steps = 25, nev = 3, max_steps = 2000)
br = continuation(prob, PALC(tangent = Bordered()), opts_br; verbosity = 0,
    bothside = false, normC = norminf)

# plot(br, plotfold=true)
####################################################################################################
# plotting function
function plotPO(x, p; k...)
    xtt = BK.get_periodic_orbit(p.prob, x, @set par_freire.β = p.p)
    plot!(xtt.t, xtt[1,:]; markersize = 2, k...)
    plot!(xtt.t, xtt[2,:]; k...)
    plot!(xtt.t, xtt[3,:]; marker = :d, markersize = 1, legend = false, k...)
end

# record function
function recordPO(x, p; k...)
    xtt = BK.get_periodic_orbit(p.prob, x, @set par_freire.β = p.p)
    period = BK.getperiod(p.prob, x, p.p)
    return (period = period, max = maximum(xtt[1,:]), min = minimum(xtt[1,:]))
end

# newton parameters
optn_po = NewtonPar(verbose = false, tol = 1e-8,  max_iterations = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.025, ds= -0.001, dsmin = 1e-4, p_max = 1.8, p_min=-5., max_steps = 130, newton_options = (@set optn_po.tol = 1e-8), detect_bifurcation = 0, plot_every_step = 3, save_sol_every_step=1,)

br_coll = continuation(
    br, 4, opts_po_cont,
    PeriodicOrbitOCollProblem(50, 4; meshadapt = true, update_section_every_step = 2);
    ampfactor = 1., δp = 0.001,
    # verbosity = 2,    plot = true,
    alg = PALC(tangent = Bordered()),
    record_from_solution = recordPO,
    plot_solution = (x, p; k...) -> begin
        plotPO(x, p; k...)
        ## add previous branch
        # plot!(br, subplot=1, putbifptlegend = false)
        end,
    finalise_solution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
        # limit the period
            return z.u[end] < 200
            true
        end,
    normC = norminf)
####################################################################################################
# homoclinic
probhom, solh = generate_hom_problem(
    setproperties(br_coll.prob.prob, meshadapt=true, K = 100),
    br_coll.sol[end].x.sol,
    BK.setparam(br_coll, br_coll.sol[end].p),
    BK.getlens(br_coll);
    update_every_step = 4,
    ϵ0 = 1e-9,
    ϵ1 = 1e-6,
    # freeparams = ((@optic _.T), (@optic _.ϵ1),)
    # freeparams = ((@optic _.T), (@optic _.ϵ0)),
    freeparams = ((@optic _.ϵ0), (@optic _.ϵ1)),
    # freeparams = ((@optic _.T),),
    )

show(probhom)
HclinicBifurcationKit.generate_homoclinic_solution(probhom.bvp, t->ones(3), 1.)
####################################################################################################
# same with shooting
using OrdinaryDiffEq
probsh = ODEProblem(freire!, copy(z0), (0., 1000.), par_freire; abstol = 1e-12, reltol = 1e-10)

# newton parameters
optn_po = NewtonPar(verbose = false, tol = 1e-8,  max_iterations = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.15, ds= -0.0001, dsmin = 1e-4, p_max = 1.8, p_min=-5., max_steps = 130, newton_options = (@set optn_po.tol = 1e-8), tol_stability = 1e-4, detect_bifurcation = 0, plot_every_step = 20, save_sol_every_step=1)

br_sh = continuation(
    br, 4, opts_po_cont,
    ShootingProblem(10, probsh, Rodas5P(); parallel = false);
    ampfactor = 1.0, δp = 0.001,
    # verbosity = 2,    plot = true,
    record_from_solution = recordPO,
    plot_solution = (x, p; k...) -> begin
        plotPO(x, p; k...)
        ## add previous branch
        # plot!(br, subplot=1, putbifptlegend = false)
        end,
    finalise_solution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
        # limit the period
            return z.u[end] < 300
            true
        end,
    normC = norminf)

# _sol = get_periodic_orbit(br_sh.prob.prob, br_sh.sol[end].x, BK.setparam(br_sh,  br_sh.sol[end].p))
        # plot(_sol.t, _sol[:,:]')
#######################################
# homoclinic
probhom, solh = generate_hom_problem(
    br_sh.prob.prob, br_sh.sol[end].x,
    BK.setparam(br_sh, br_sh.sol[end].p),
    BK.getlens(br_sh);
    verbose = true,
    update_every_step = 4,
    ϵ0 = 7e-8,
    ϵ1 = 8e-8,
    # freeparams = ((@optic _.T), (@optic _.ϵ1),)
    # freeparams = ((@optic _.T), (@optic _.ϵ0)),
    freeparams = ((@optic _.ϵ0), (@optic _.ϵ1)), # WORK BEST
    # freeparams = ((@optic _.T),),
    )

solh.x[2] .=0

show(probhom)