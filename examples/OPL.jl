using Revise, Plots
using LinearAlgebra, Test, ForwardDiff
using BifurcationKit, Test
using HclinicBifurcationKit
const BK = BifurcationKit

recordFromSolution(x, p;k...) = (D₂₃ = x[6], β = x[1],)
####################################################################################################
function OPL!(dz, u, p, t)
    (;b, σ, g, a, D₂₁⁰, D₂₃⁰)  = p
    β, p₂₁, p₂₃, p₃₁, D₂₁, D₂₃ = u
    dz[1] = -σ * β + g * p₂₃
    dz[2] = -p₂₁ - β * p₃₁ + a * D₂₁
    dz[3] = -p₂₃ + β * D₂₃ - a * p₃₁
    dz[4] = -p₃₁ + β * p₂₁ + a * p₂₃
    dz[5] = -b * (D₂₁ - D₂₁⁰) - 4a * p₂₁ - 2β * p₂₃
    dz[6] = -b * (D₂₃ - D₂₃⁰) - 2a * p₂₁ - 4β * p₂₃
    dz
end

OPL(z, p) = OPL!(similar(z), z, p, 0)
par_OPL = (b = 1.2, σ = 2.0, g=50., a = 1., D₂₁⁰ = -1., D₂₃⁰ = 0.)
z0 = zeros(6)
prob = BK.BifurcationProblem(OPL, z0, par_OPL, (@optic _.a); record_from_solution = recordFromSolution)

opts_br = ContinuationPar(p_min = -1., p_max = 8., ds = 0.001, dsmax = 0.06, n_inversion = 6, detect_bifurcation = 3, max_bisection_steps = 25, nev = 6, plot_every_step = 20, max_steps = 100, save_sol_every_step = 1, detect_fold = true)
br = continuation(prob, PALC(tangent = Secant()), opts_br;
    bothside = false, normC = norminf)

plot(br, plotfold=true)

br2 = continuation(br, 1)

plot(br, br2)
####################################################################################################
sn_br = continuation(br, 1, (@optic _.b), ContinuationPar(opts_br, detect_bifurcation = 1, save_sol_every_step = 1, max_steps = 80) ;
    alg = PALC(),
    verbosity = 0,
    detect_codim2_bifurcation = 2,
    start_with_eigen = true,
    update_minaug_every_step = 1,
    bothside = true,
    )

plot(sn_br)
bt = get_normal_form(sn_br, 2, verbose = true, detailed = true, autodiff = false)

hopf_br = continuation(br, 2, (@optic _.b), ContinuationPar(opts_br, detect_bifurcation = 1, save_sol_every_step = 1, max_steps = 140),
    detect_codim2_bifurcation = 2,
    start_with_eigen = true,
    update_minaug_every_step = 1,
    bothside = true,
    )

hopf_br2 = continuation(br2, 1, (@optic _.b), ContinuationPar(opts_br, detect_bifurcation = 1, save_sol_every_step = 1, max_steps = 140),
    detect_codim2_bifurcation = 2,
    start_with_eigen = true,
    update_minaug_every_step = 1,
    bothside = true,
    )

plot(sn_br, vars = (:a, :b), branchlabel = ["SN"], )
plot!(hopf_br, branchlabel = ["Hopf"], vars = (:a, :b))
plot!(hopf_br2, branchlabel = ["Hopf"], vars = (:a, :b))
ylims!(0,1.5)

####################################################################################################
# plotting function
function plotPO(x, p; k...)
    xtt = BK.get_periodic_orbit(p.prob, x, set(getparams(p.prob), BK.getlens(p.prob), p.p))
    plot!(xtt.t, xtt[1,:]; markersize = 2, k...)
    plot!(xtt.t, xtt[6,:]; k...)
    scatter!(xtt.t, xtt[1,:]; markersize = 1, legend = false, k...)
end

# record function
function recordPO(x, p; k...)
    xtt = BK.get_periodic_orbit(p.prob, x, set(getparams(p.prob), BK.getlens(p.prob), p.p))
    period = BK.getperiod(p.prob, x, p.p)
    return (max = maximum(xtt[6,:]), min = minimum(xtt[6,:]), period = period, )
end

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  max_iterations = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.05, ds= 0.001, dsmin = 1e-4, p_max = 6.8, p_min=-5., max_steps = 100, newton_options = (@set optn_po.tol = 1e-8), detect_bifurcation = 0, plot_every_step = 3, save_sol_every_step=1,)

br_coll = continuation(
    # br, 2,
    br2, 1,
    opts_po_cont,
    PeriodicOrbitOCollProblem(40, 4; meshadapt = true, update_section_every_step = 2);
    ampfactor = 1., δp = 0.0015,
    verbosity = 2,    plot = true,
    alg = PALC(tangent = Bordered()),
    record_from_solution = recordPO,
    plot_solution = (x, p; k...) -> begin
        plotPO(x, p; k...)
        ## add previous branch
        plot!(br, subplot=1, putbifptlegend = false)
        plot!(br2, subplot=1, putbifptlegend = false)
        end,
    finalise_solution = (z, tau, step, contResult; prob = nothing, kwargs...) -> begin
        # limit the period
            return z.u[end] < 150
            true
        end,
    normC = norminf)

_sol = get_periodic_orbit(br_coll, length(br_coll))
BK.plot(_sol.t, _sol.u'; marker = :d, markersize = 1, title = "Last periodic orbit on branch")
####################################################################################################
# homoclinic
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
    # freeparams = ((@optic _.T), (@optic _.ϵ1),)
    # freeparams = ((@optic _.T), (@optic _.ϵ0)),
    freeparams = ((@optic _.ϵ0), (@optic _.ϵ1)), # WORK BEST
    # freeparams = ((@optic _.T),),
    testOrbitFlip = false,
    testInclinationFlip = false
    )

#####

_sol = get_homoclinic_orbit(probhom, solh, BK.getparams(probhom);)
plot(_sol.t, _sol[:,:]', marker = :d, markersize = 1, title = "Initial guess for homoclinic orbit")

optn_hom = NewtonPar(verbose = true, tol = 1e-10, max_iterations = 5)
optc_hom = ContinuationPar(newton_options = optn_hom, ds = -0.0001, dsmin = 1e-5, plot_every_step = 10, max_steps = 100, detect_bifurcation = 0, detect_event = 2, save_sol_every_step = 1, p_min = -1.01)

br_hom_c = continuation(
            deepcopy(probhom), solh, (@optic _.b),
            PALC(tangent = Bordered()),
            # PALC(),
            # MoorePenrose(),
            setproperties(optc_hom, max_steps = 100, save_sol_every_step = 1, dsmax = 4e-2, plot_every_step = 10, p_max = 1.5);
    verbosity = 1, plot = true,
    # callback_newton = BK.cbMaxNorm(1e1),
    # bothside = true,
    normC = norminf,
    plot_solution = (x,p;k...) -> begin
        𝐇𝐨𝐦 = p.prob
        par0 = set(BK.getparams(𝐇𝐨𝐦), BK.getlens(𝐇𝐨𝐦), x.x[end][1])
        par0 = set(par0, (@optic _.b), p.p)
        sol = get_homoclinic_orbit(𝐇𝐨𝐦, x, par0)
        m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
        plot!(sol.t, sol[:,:]',subplot=3, markersize = 1, marker=m)
    end,
    )


br_hom_c.branch |> vscodedisplay

plot(sn_br, vars = (:a, :b), branchlabel = "SN", )
plot!(hopf_br, branchlabel = "AH₀", vars = (:a, :b))
plot!(hopf_br2, branchlabel = "AH₁₂", vars = (:a, :b))
plot!(br_hom_c, branchlabel = "H₀", vars = (:a, :b))
ylims!(0,1.5)
####################################################################################################
# same with shooting
using OrdinaryDiffEq
probsh = ODEProblem(OPL!, copy(z0), (0., 1000.), par_OPL; abstol = 1e-12, reltol = 1e-10)

# newton parameters
optn_po = NewtonPar(verbose = true, tol = 1e-8,  max_iterations = 25)

# continuation parameters
opts_po_cont = ContinuationPar(dsmax = 0.075, ds= -0.001, dsmin = 1e-4, p_max = 6.8, p_min=-5., max_steps = 130, newton_options = (@set optn_po.tol = 1e-8), tol_stability = 1e-4, detect_bifurcation = 0, plot_every_step = 10, save_sol_every_step=1)

br_sh = continuation(
    # br, 2,
    br2, 1,
    opts_po_cont,
    ShootingProblem(8, probsh, Rodas5P(); parallel = true, abstol = 1e-13, reltol = 1e-11);
    ampfactor = 1., δp = 0.0015,
    verbosity = 2,    plot = true,
    record_from_solution = recordPO,
    # alg = MoorePenrose(),
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

_sol = get_periodic_orbit(br_sh, length(br_sh))
plot(_sol)

#######################################
# homoclinic
probhom, solh = generate_hom_problem(
    br_sh.prob.prob, br_sh.sol[end].x,
    BK.setparam(br_sh, br_sh.sol[end].p),
    BK.getlens(br_sh);
    verbose = true,
    update_every_step = 4,
    # ϵ0 = 1e-6, ϵ1 = 1e-5,
    t0 = 75, t1 = 25,
    # freeparams = ((@optic _.T), (@optic _.ϵ1),)
    # freeparams = ((@optic _.T), (@optic _.ϵ0)),
    # freeparams = ((@optic _.ϵ0), (@optic _.ϵ1)), # WORK BEST
    freeparams = ((@optic _.T),),
    )

_sol = get_homoclinic_orbit(probhom, solh, BK.getparams(probhom); saveat=.1)
plot(plot(_sol[1,:], _sol[2,:]), plot(_sol.t, _sol[1:4,:]'))

optn_hom = NewtonPar(verbose = true, tol = 1e-9, max_iterations = 7)
    optc_hom = ContinuationPar(newton_options = optn_hom, ds = -1e-4, dsmin = 1e-6, dsmax = 1e-3, plot_every_step = 1, max_steps = 10, detect_bifurcation = 0, save_sol_every_step = 1)

br_hom_sh = continuation(
            deepcopy(probhom), solh, (@optic _.b),
            # PALC(tangent = Bordered()),
            PALC(),
            # ANM(6, 1e-8)
            # MoorePenrose(),
            setproperties(optc_hom, max_steps = 600, save_sol_every_step = 1, dsmax = 12e-2, plot_every_step = 3, p_max = 7., detect_event = 2, a = 0.9);
    verbosity = 3, plot = true,
    callback_newton = BK.cbMaxNorm(1e0),
    normC = norminf,
    plot_solution = (x,p;k...) -> begin
        𝐇𝐨𝐦 = p.prob
        par0 = set(BK.getparams(𝐇𝐨𝐦), BK.getlens(𝐇𝐨𝐦), x.x[end][1])
        par0 = set(par0, (@optic _.b), p.p)
        sol = get_homoclinic_orbit(𝐇𝐨𝐦, x, par0)
        m = (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem && 𝐇𝐨𝐦.bvp.meshadapt) ? :d : :none
        plot!(sol.t, sol[1:6,:]',subplot=3, markersize = 1, marker=m)
    end,
    )

_sol = get_homoclinic_orbit(probhom, br_hom_sh.sol[end].x, BK.setparam(br_hom_sh, br_hom_sh.sol[end].p); saveat=.1)
plot(_sol.t, _sol[1:6,:]')


plot(sn_br, vars = (:a, :b), branchlabel = "SN", )
plot!(hopf_br, branchlabel = "AH₀", vars = (:a, :b))
plot!(hopf_br2, branchlabel = "AH₁₂", vars = (:a, :b))
plot!(br_hom_c, branchlabel = "Hc₀", vars = (:a, :b))
plot!(br_hom_sh, branchlabel = "Hsh₀", vars = (:a, :b), linewidth = 3)
ylims!(0,1.5)
