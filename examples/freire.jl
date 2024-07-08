using Revise, Plots
using LinearAlgebra, Test
using BifurcationKit, Test
using HclinicBifurcationKit
const BK = BifurcationKit

recordFromSolution(x, p) = (x = x[1], y = x[2])
####################################################################################################
function freire!(dz, u, p, t)
    (;ν, β, A₃, B₃, r, ϵ) = p
    x, y, z = u
    dz[1] = (-ν*x + β*(y-x) - A₃*x^3 + B₃*(y-x)^3 + ϵ)/r
    dz[2] =    -β*(y-x) - z - B₃*(y-x)^3
    dz[3] = y
    dz
end

freire(z, p) = freire!(similar(z), z, p, 0)
par_freire = (ν = -0.75, β = -0.1, A₃ = 0.328578, B₃ = 0.933578, r = 0.6, ϵ = 0.01)
z0 = [0.7,0.3,0.1]
z0 = zeros(3)
prob = BK.BifurcationProblem(freire, z0, par_freire, (@lens _.β); record_from_solution = recordFromSolution)

opts_br = ContinuationPar(p_min = -1.4, p_max = 2.8, ds = 0.001, dsmax = 0.05, n_inversion = 6, max_bisection_steps = 25, nev = 3, max_steps = 2000)
br = continuation(prob, PALC(tangent = Bordered()), opts_br;
bothside = false, normC = norminf)

plot(br, plotfold=true)
####################################################################################################
# DefCon
alg = DefCont(deflation_operator = DeflationOperator(2, dot, .001, [z0]), perturb_solution = (x,p,id) -> (x  .+ 0.03 .* rand(length(x))))

opts_br_dc = ContinuationPar(p_min = -1.4, p_max = 2.8, ds = 0.001, dsmax = 0.01, detect_bifurcation = 3, max_bisection_steps = 25, nev = 3, plot_every_step = 100, max_steps = 1000)
@set! opts_br_dc.newton_options.verbose = false
brdc = continuation(prob, alg, opts_br_dc; verbosity = 0, normC = norminf)
####################################################################################################
sn_br = continuation(br, 2, (@lens _.ν), ContinuationPar(opts_br, detect_bifurcation = 1, save_sol_every_step = 1, dsmax = 0.01, max_steps = 80) ;
    alg = PALC(),
    detect_codim2_bifurcation = 2,
    start_with_eigen = true,
    update_minaug_every_step = 1,
    bothside = true,
    )

plot(sn_br)
bt = get_normal_form(sn_br, 2, verbose = true, detailed = true, autodiff = false)

hopf_br = continuation(br, 4, (@lens _.ν), ContinuationPar(opts_br, detect_bifurcation = 1, save_sol_every_step = 1, max_steps = 140, dsmax = 0.02, n_inversion = 6),
    detect_codim2_bifurcation = 2,
    start_with_eigen = true,
    update_minaug_every_step = 1,
    bothside = true,
    )

plot(sn_br, hopf_br, ylims = (-0.1, 1.25))
####################################################################################################
function plotHom(x,p;k...)
    𝐇𝐨𝐦 = p.prob
    par0 = set(BK.getparams(𝐇𝐨𝐦), BK.getlens(𝐇𝐨𝐦), x.x[end][1])
    par0 = set(par0, p.lens, p.p)
    sol = get_homoclinic_orbit(𝐇𝐨𝐦, x, par0)
    m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
    plot!(sol.t, sol[1:3,:]',subplot=3, markersize = 1, marker=m)
end

btpt = get_normal_form(sn_br, 2; nev = 3, autodiff = false)

br_hom_c = continuation(
            prob,
            btpt,
            PeriodicOrbitOCollProblem(50, 3; meshadapt = false, K = 200),
            PALC(tangent = Bordered()),
            ContinuationPar(opts_br, max_steps = 30, save_sol_every_step = 1, dsmax = 1e-2, plot_every_step = 1, p_min = -1.01, ds = 0.001, detect_event = 2, detect_bifurcation = 0);
    verbosity = 1, plot = false,
    ϵ0 = 1e-5, amplitude = 2e-3,
    # freeparams = ((@lens _.T), (@lens _.ϵ1),)
    # freeparams = ((@lens _.T), (@lens _.ϵ0)),
    freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)),
    normC = norminf,
    plot_solution = plotHom,
    update_every_step = 4,
    )

plot(sn_br, hopf_br, ylims = (0, 1.25), branchlabel = ["SN", "HOPF"])
plot!(br_hom_c, branchlabel = "Hom")
####################################################################################################
# plotting function
function plotPO(x, p; k...)
    xtt = BK.get_periodic_orbit(p.prob, x, p.p)
    plot!(xtt.t, xtt[1,:]; markersize = 2, k...)
    plot!(xtt.t, xtt[2,:]; k...)
    plot!(xtt.t, xtt[3,:]; marker = :d, markersize = 1, legend = false, k...)
end

# record function
function recordPO(x, p)
    xtt = BK.get_periodic_orbit(p.prob, x, p.p)
    period = BK.getperiod(p.prob, x, p.p)
    return (period = period, max = maximum(xtt[1,:]), min = minimum(xtt[1,:]))
end

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  max_iterations = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, p_max = 1.8, p_min=-5., max_steps = 130, newton_options = (@set optn_po.tol = 1e-8), detect_bifurcation = 0, plot_every_step = 3)

br_coll = continuation(
    br, 4, opts_po_cont,
    PeriodicOrbitOCollProblem(50, 4; meshadapt = false, update_section_every_step = 2);
    ampfactor = 1., δp = 0.001,
    verbosity = 2,    plot = true,
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

_sol = get_periodic_orbit(br_coll, length(br_coll))
plot(_sol.t, _sol[:,:]')
####################################################################################################
# homoclinic
probhom, solh = generate_hom_problem(
    setproperties(br_coll.prob.prob, meshadapt=true, K = 100),
    br_coll.sol[end].x,
    BK.setparam(br_coll, br_coll.sol[end].p),
    BK.getlens(br_coll);
    update_every_step = 4,
    ϵ0 = 1e-9,
    ϵ1 = 1e-6,
    # freeparams = ((@lens _.T), (@lens _.ϵ1),)
    # freeparams = ((@lens _.T), (@lens _.ϵ0)),
    freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)),
    # freeparams = ((@lens _.T),),
    )

#####

_sol = get_homoclinic_orbit(probhom, solh, BK.getparams(probhom);)
plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[:,:]'))

optn_hom = NewtonPar(verbose = true, tol = 1e-10, max_iterations = 5)
optc_hom = ContinuationPar(newton_options = optn_hom, ds = 0.0001, dsmin = 1e-5, plot_every_step = 10, max_steps = 100, detect_bifurcation = 0, detect_event = 2, save_sol_every_step = 1, p_min = -1.01)

solh.x[2] .=0

br_hom_c = continuation(
            deepcopy(probhom), solh, (@lens _.ν),
            PALC(tangent = Bordered()),
            ContinuationPar(optc_hom, max_steps = 130, save_sol_every_step = 1, dsmax = 1e-2, plot_every_step = 1);
    verbosity = 4, plot = true,
    normC = norminf,
    plot_solution = plotHom,
    )

using PrettyTables
br_hom_c.branch[end-20:end] |> pretty_table

plot(sn_br, hopf_br, ylims = (0, 1.25))
plot!(br_hom_c)

_sol = get_homoclinic_orbit(probhom, br_hom_c.sol[end].x, BK.setparam(br_hom_c,  br_hom_c.sol[end].p))
plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[:,:]'))
####################################################################################################
# same with shooting
using OrdinaryDiffEq
probsh = ODEProblem(freire!, copy(z0), (0., 1000.), par_freire; abstol = 1e-12, reltol = 1e-10)

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  max_iterations = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.15, ds= -0.0001, dsmin = 1e-4, p_max = 1.8, p_min=-5., max_steps = 130, newton_options = (@set optn_po.tol = 1e-8), tol_stability = 1e-4, detect_bifurcation = 0, plot_every_step = 20, save_sol_every_step=1)

br_sh = continuation(
    br, 4, opts_po_cont,
    ShootingProblem(10, probsh, Rodas5P(); parallel = false);
    ampfactor = 1.0, δp = 0.001,
    verbosity = 2,    plot = true,
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

_sol = get_periodic_orbit(br_sh.prob.prob, br_sh.sol[end].x, BK.setparam(br_sh,  br_sh.sol[end].p))
plot(_sol.t, _sol[1:3,:]')
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
    # freeparams = ((@lens _.T), (@lens _.ϵ1),)
    # freeparams = ((@lens _.T), (@lens _.ϵ0)),
    freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)), # WORK BEST
    # freeparams = ((@lens _.T),),
    )

solh.x[2] .= 0

_sol = get_homoclinic_orbit(probhom, solh, BK.getparams(probhom))
plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[1:3,:]'))

optn_hom = NewtonPar(verbose = true, tol = 1e-10, max_iterations = 7)
optc_hom = ContinuationPar(newton_options = optn_hom, ds = 1e-4, dsmin = 1e-6, dsmax = 1e-3, plot_every_step = 1,max_steps = 10, detect_bifurcation = 0, save_sol_every_step = 1)

br_hom_sh = continuation(
            deepcopy(probhom), solh, (@lens _.ν),
            PALC(tangent = Bordered()),
            # MoorePenrose(),
            ContinuationPar(optc_hom, max_steps = 30, save_sol_every_step = 1, dsmax = 3e-2, plot_every_step = 10, detect_event = 2);
    verbosity = 4, plot = true,
    callback_newton = BK.cbMaxNorm(1e1),
    normC = norminf,
    plot_solution = plotHom,
    )

_sol = get_homoclinic_orbit(br_hom_sh.prob.VF.F, br_hom_sh.sol[end].x, BK.setparam(br_hom_sh, br_hom_sh.sol[end].p))
plot(_sol.t, _sol[1:3,:]')

plot(hopf_br, br_hom_c, br_hom_sh)

#######
br_hom_sh = continuation(
            prob,
            btpt,
            ShootingProblem(12, probsh, Rodas5P(); parallel = true, abstol = 1e-13, reltol = 1e-12),
            PALC(tangent = Bordered()),
            ContinuationPar(optc_hom, max_steps = 15000, save_sol_every_step = 1, ds = 1e-3, dsmax = 3e-2, plot_every_step = 50, detect_event = 2, a = 0.9, p_min = -1.01);
    verbosity = 1, plot = true,
    ϵ0 = 1e-6, amplitude = 2e-2,
    update_every_step = 2,
    # freeparams = ((@lens _.T), (@lens _.ϵ1),),
    # freeparams = ((@lens _.T), (@lens _.ϵ0)),
    # freeparams = ((@lens _.ϵ0), (@lens _.ϵ1)),
    freeparams = ((@lens _.T),),
    normC = norminf,
    plot_solution = plotHom,
    maxT = 45.,
    callback_newton = BK.cbMaxNorm(1e1),
    )
