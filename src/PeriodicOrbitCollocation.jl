getF(hom::HomoclinicHyperbolicProblemPBC{Tbvp}, x, p) where {Tbvp <: PeriodicOrbitOCollProblem} = BK.residual(hom.bvp.prob_vf, x, p)

get_homoclinic_orbit(hom::HomoclinicHyperbolicProblemPBC{Tbvp}, x::ArrayPartition, par; k...) where {Tbvp <: PeriodicOrbitOCollProblem} = get_periodic_orbit(hom.bvp, vcat(x.x[1], hom.T), par)

function generate_hom_solution(pb::PeriodicOrbitOCollProblem, orbit0, T)
    orbit = t -> orbit0(-T + t * (2T))
    generate_solution(pb, orbit, 1.)[1:end-1]
end

function initBVPforPBC(bvp::PeriodicOrbitOCollProblem, prob_vf, Hom; N, T, ϵ)
    @reset bvp.N = N
    bvp = setproperties(bvp; prob_vf = prob_vf, ϕ = zeros(length(bvp)), xπ = zeros(length(bvp)), update_section_every_step = 0)
    _N, m, Ntst = size(bvp)
    bvp = BK.set_collocation_size(bvp, Ntst, m)
    cache = BK.POCollCache(eltype(bvp), Ntst, N, m)
    @reset bvp.cache = cache
    xflow = generate_hom_solution(bvp, t -> Hom.orbit(t, ϵ), T)
    BK.updatesection!(bvp, vcat(xflow, 2T), BK.getparams(bvp))
    return xflow, bvp
end

"""
$(SIGNATURES)

This function generates an initial guess for the solution of the problem `pb` based on the orbit `t -> orbit(t)` for t ∈ [-T,T] and half time return `T`.
"""
function generate_homoclinic_solution(pb::PeriodicOrbitOCollProblem, orbit, T)
    n, _m, Ntst = size(pb)
    ts = BK.get_times(pb)
    Nt = length(ts)
    ci = zeros(eltype(pb), n, Nt)
    for (l, t) in pairs(ts)
        ci[:, l] .= orbit(-T + t * (2T))
    end
    return vec(ci)
end

"""
Implements
    ∫ < u - v, vₜ >
"""
@views function phase_condition_PBC(pb::PeriodicOrbitOCollProblem, (u, uc), (L, ∂L))
    Ty = eltype(uc)
    phase = zero(Ty)

    n, m, Ntst = size(pb)

    guj = zeros(Ty, n, m)
    uj  = zeros(Ty, n, m+1)

    vc = BK.get_time_slices(pb.ϕ, size(pb)...)
    gvj = zeros(Ty, n, m)
    gdvj = zeros(Ty, n, m)
    vj  = zeros(Ty, n, m+1)

    ω = pb.mesh_cache.gauss_weight

    rg = UnitRange(1, m+1)
    @inbounds for j in 1:Ntst
        uj .= uc[:, rg]
        vj .= vc[:, rg]
        mul!(guj, uj, L)
        mul!(gvj, vj, L)
        mul!(gdvj, vj, ∂L)
        @inbounds for l in 1:m
            phase += dot(guj[:, l], gdvj[:, l]) * ω[l]
            phase -= dot(gvj[:, l], gdvj[:, l]) * ω[l]
        end
        rg = rg .+ m
    end
    return phase / getperiod(pb, u, nothing)
end

@views function (hom::HomoclinicHyperbolicProblemPBC{Tbvp, Nf})(x::ArrayPartition, par0) where {Tbvp <: PeriodicOrbitOCollProblem, Nf}
    @unpack N = hom
    coll = hom.bvp
    ns = hom.nStable
    nu = hom.nUnstable

    _u = x.x[1]         # orbit
    xsaddle = x.x[2]    # saddle point
    Ys = x.x[3]         #   stable part for CIS algo
    Yu = x.x[4]         # unstable part for CIS algo
    # get homoclinic parameters
    T, ϵ0, ϵ1 = _changeHomParameters(hom, x.x[5])

    @assert size(Ys) == (N - ns, ns) "size(Ys) = $(size(Ys)) != $((N - ns, ns))"
    @assert size(Yu) == (N - nu, nu) "size(Yu) = $(size(Yu)) != $((N - nu, nu))"

    # get the updated parameter
    param = x.x[5][1]
    lens = hom.lens
    par = set(par0, lens, param)

    # we hack the functional for periodic orbits
    u = vcat(_u, T)
    uc = BK.get_time_slices(coll, u)
    x0 = uc[:, 1]
    x1 = uc[:, end]

    # version of collocation problem without boundary condition
    _resuc = similar(uc, N, size(uc,2)-1)
    resu = vec(_resuc)
    BK.functional_coll_bare!(coll, _resuc, uc, T, BK.get_Ls(coll.mesh_cache), par)

    # F(xsaddle, par) = 0
    Fx = getF(hom, xsaddle, par)

    # ricatti equations
    J = ForwardDiff.jacobian(x -> getF(hom, x, par), xsaddle)

    Tb  = ricattiBlocks(hom.Qu0, J, hom.nUnstable)
    riU = ricattiEq(Tb, Yu)

    Tb  = ricattiBlocks(hom.Qs0, J, hom.nStable)
    riS = ricattiEq(Tb, Ys)

    # projector on stable / unstable manifold
    Qu1⊥ = hom.Qu0 * vcat(-Yu', I(size(Yu,1)))
    uP = Qu1⊥' * (x0 - xsaddle)
    Qs1⊥ = hom.Qs0 * vcat(-Ys', I(size(Ys,1)))
    sP = Qs1⊥' * (x1 - xsaddle)

    # set distance to saddle
    outnrm = zeros(eltype(x0), 1 + Nf)
    outnrm[1] = norm(x0 .- xsaddle) - ϵ0
    outnrm[2] = norm(x1 .- xsaddle) - ϵ1
    if Nf == 2
        outnrm[3] = phase_condition_PBC(coll, (u, uc), BK.get_Ls(coll.mesh_cache))
    end

    out = ArrayPartition(resu, Fx, riU, riS, uP, sP, outnrm)
    return out
end

using SciMLBase: AbstractTimeseriesSolution
"""
$(SIGNATURES)

Generate a homoclinic to hyperbolic saddle problem from a periodic solution obtained with problem `pb`.

## Arguments
- `coll` a `PeriodicOrbitOCollProblem` which provide basic information, like the number of time slices `M`
- `x::AbstractArray` initial guess
- `pars` parameters
- `lensHom::Lens` parameter axis for continuation
- `ϵ0, ϵ1`: specify the distance to the saddle point of x₀, x₁
- `t0, t1`: specify the time corresponding to x₀, x₁. Overwrite the part with `ϵ0, ϵ1` if set.

## Optional arguments
You can pass the same arguments to the constructor of `::HomoclinicHyperbolicProblemPBC`.

## Output
- returns a `HomoclinicHyperbolicProblemPBC` and an initial guess.
"""
function generate_hom_problem(coll::PeriodicOrbitOCollProblem,
                            x::AbstractArray,
                            pars,
                            lensHom::BK.AllOpticTypes;
                            verbose = false,
                            ϵ0 = 1e-5, ϵ1 = 1e-5,
                            t0 = 0, t1 = 0,
                            maxT = Inf,
                            freeparams = ((@optic _.ϵ0), (@optic _.T)),
                            kw...)
    println("="^40)
    @assert coll.N > 0
    T = getperiod(coll, x)
    time = BK.get_times(coll) .* T
    xc = BK.get_time_slices(coll, x)
    indmax = size(xc, 2)

    # convert solution to homogenous mesh
    solpo = BK.POSolution(deepcopy(coll), x)

    # find the saddle point as minimum of vector field norm
    xc = BK.get_time_slices(coll, x)
    ind_saddle = argmin(norm(BK.residual(coll.prob_vf, xc[:, i], pars)) for i = 1:indmax)
    xsaddle = xc[:, ind_saddle]
    tsaddle = time[ind_saddle]
    BK._newton(coll.prob_vf, xsaddle, pars, NewtonPar(verbose = true))

    if t1 == t0 == 0
        # find x0 and x1 on the unstable / stable subspace
        indUS = findfirst(norm(solpo(t) - xsaddle) > ϵ0 for t in time .+ tsaddle)
        t0 = mod(time[indUS]+tsaddle, T)
        x0 = solpo(t0)
        indS = findlast(norm(solpo(t) - xsaddle) > ϵ1 for t in time .+ t0)
        t1 = time[indS]+t0
        x1 = solpo(t1)
    else
        x0 = solpo(t0)
        x1 = solpo(t1)
        indUS, indS = 0, 0
    end

    # we put a uniform mesh in bvp even if coll is non uniform
    n, m, Ntst = size(coll)
    bvp = deepcopy(coll)
    bvp = BK.set_collocation_size(bvp, Ntst, m)
    @reset bvp.update_section_every_step = 0
    BK.update_mesh!(bvp, LinRange{eltype(coll)}( 0, 1, Ntst + 1) |> collect)
    bvp = BK._set_params_po(bvp, pars)

    Thom = min(mod(t1-t0, T), maxT)
    xflow = mapreduce(t -> solpo(t0 + t * Thom), vcat, BK.get_times(bvp))
    BK.updatesection!(bvp, vcat(xflow, Thom), BK.getparams(bvp))

    # create Homoclinic parameters
    ϵ0hom = norm(x0 - xsaddle)
    ϵ1hom = norm(x1 - xsaddle)

    # define problem for Homoclinic functional
    J = BK.jacobian(coll.prob_vf, xsaddle, pars)
    𝐇𝐨𝐦 = HomoclinicHyperbolicProblemPBC(bvp, lensHom, length(xsaddle), copy(J);  ϵ0 = ϵ0hom, ϵ1 = ϵ1hom, T = Thom, freeparams = freeparams, kw...)

    @assert BK.getparams(𝐇𝐨𝐦) == pars "Errors with setting the parameters. Please an issue on the website of BifurcationKit."

    if verbose
        println("┌─ tsaddle  = $tsaddle")
        println("├─ t0       = $t0")
        println("├─ t1       = $t1")
        println("├─ T        = $Thom")
        println("└─ is,i0,i1 = $((ind_saddle, indUS, indS))")
    end

    ns = 𝐇𝐨𝐦.nStable
    nu = 𝐇𝐨𝐦.nUnstable
    p1 = BK._get(pars, lensHom)

    xhom = ArrayPartition(xflow,
        xsaddle,
        zeros(eltype(xsaddle), n - ns, ns),
        zeros(eltype(xsaddle), n - nu, nu),
        [p1, map(x -> BK._get(𝐇𝐨𝐦,x), freeparams)...]
        )

    return 𝐇𝐨𝐦, xhom, pars, xhom
end

function generate_hom_problem(coll::PeriodicOrbitOCollProblem,
                            x::NamedTuple{(:mesh, :sol, :_mesh), Tuple{Vector{Tp}, Vector{Tp}, Vector{Tp}}},
                            pars,
                            lensHom::BK.AllOpticTypes; k...) where Tp
    n, m, _ = size(coll)
    coll2 = deepcopy(coll)
    BK.update_mesh!(coll2, x.mesh[1:m:end])
    generate_hom_problem(coll2, x.sol, pars, lensHom; k...)
end
