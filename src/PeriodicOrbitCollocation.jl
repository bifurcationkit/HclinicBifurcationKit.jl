getF(hom::HomoclinicHyperbolicProblemPBC{Tdisc}, x, p) where {Tdisc <: Collocation} = BK.residual(BK.get_discretization(hom).prob_vf, x, p)

"""
$(TYPEDSIGNATURES)

Reconstruct the homoclinic orbit from a solution `x` of the homoclinic problem `hom`, for example a point stored on a branch of homoclinic orbits (`br.sol[ind].x`). `par` are the parameters at which the solution was computed (e.g. `BK.setparam(br, br.sol[ind].p)`).

The extra `kwargs` (e.g. `saveat`) are forwarded to the ODE solver used to integrate the orbit.

# Output
Returns the homoclinic orbit as a time solution with fields `t` (time mesh) and `u` (states), which can be plotted or interpolated in time.
"""
get_homoclinic_orbit(hom::HomoclinicHyperbolicProblemPBC{Tdisc}, x::ArrayPartition, par; k...) where {Tdisc <: Collocation} = get_periodic_orbit(BK.get_discretization(hom), vcat(x.x[1], hom.T), par)

function generate_hom_solution(pb::Collocation, orbit0, T)
    orbit = t -> orbit0(-T + t * (2T))
    generate_solution(pb, orbit, 1.)[1:end-1]
end

function initBVPforPBC(coll::Collocation, prob_vf, Hom; N, T, ϵ)
    @reset coll.N = N
    coll = setproperties(coll; prob_vf = prob_vf, ϕ = zeros(length(coll)), xπ = zeros(length(coll)), update_section_every_step = 0)
    _N, m, Ntst = size(coll)
    coll = BK.set_collocation_size(coll, Ntst, m)
    cache = BK.POCollCache(eltype(coll), Ntst, N, m)
    @reset coll.cache = cache
    xflow = generate_hom_solution(coll, t -> Hom.orbit(t, ϵ), T)
    BK.updatesection!(coll, vcat(xflow, 2T), BK.getparams(coll))
    return xflow, coll
end

"""
$(SIGNATURES)

This function generates an initial guess for the solution of the problem `pb` based on the orbit `t -> orbit(t)` for t ∈ [-T, T] and half time return `T`.
"""
function generate_homoclinic_solution(disc::Collocation, orbit, T)
    n, _m, Ntst = size(disc)
    ts = BK.get_times(disc)
    Nt = length(ts)
    ci = zeros(eltype(disc), n, Nt)
    for (l, t) in pairs(ts)
        ci[:, l] .= orbit(-T + t * (2T))
    end
    return vec(ci)
end

"""
Implements
    ∫ < u - v, vₜ >
"""
@views function phase_condition_PBC(pb::Collocation, (u, uc), (L, ∂L))
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

# residual function
@views function (hom::HomoclinicHyperbolicProblemPBC{Tdisc, Nf})(x::ArrayPartition, par0) where {Tdisc <: Collocation, Nf}
    (; N) = hom
    coll = BK.get_discretization(hom)
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
    _resuc = similar(uc, N, size(uc, 2) - 1)
    resu = vec(_resuc)
    BK.po_residual_bare!(coll, _resuc, uc, T, BK.get_Ls(coll.mesh_cache), par)

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
    _select_homoclinic_arc(τs, ft, m, Ntst, isaddle, t0, t1, T)

Select, on the stored points of a periodic-orbit solution (no interpolation), the
contiguous block of whole collocation intervals used as an initial guess for the truncated
homoclinic orbit.

- `τs`: interval boundaries (length `Ntst + 1`)
- `ft`: node fractions of the periodic orbit (length `1 + m * Ntst`)
- `isaddle`: index of the node closest to the saddle point

If `t0 == t1 == 0`, the single interval containing the closest node is removed.
Otherwise the whole intervals covering the window `[t0, t1]` are kept.

Return a NamedTuple with:
- `kept`: indices (in `1:Ntst`) of the intervals kept along the arc,
- `uidx`: indices of the kept time slices, in arc order,
- `ulifts`: unwrapped fractions of these points,
- `widths`: widths of the kept intervals,
- `a0`: physical fraction of the arc start (a boundary),
- `L`: arc length, as a fraction of the period.
"""
function _select_homoclinic_arc(τs, ft, m, Ntst, isaddle, t0, t1, T)
    # ━━━━━━ choose the contiguous block of whole intervals ━━━━━━
    if t1 == t0 == 0
        # remove the interval containing the closest node
        if Ntst == 1
            g = 1
            kept = [1]
        else
            jrem = isaddle == 1 ? 1 : max(1, cld(isaddle - 1, m))
            g = mod1(jrem + 1, Ntst)                    # first interval after the removed one
            kept = [mod1(jrem + kk, Ntst) for kk in 1:(Ntst - 1)]
        end
    else
        f0 = mod(t0, T) / T
        f1 = mod(t1, T) / T
        # first interval of the arc: down-snap of f0 to a boundary
        g = clamp(searchsortedlast(τs, f0), 1, Ntst)
        # last interval of the arc: up-snap (interval ending at, or containing, f1)
        jlast = if f1 >= 1 - 1e-12 || f1 < 1e-12
            Ntst
        else
            kk = searchsortedfirst(τs, f1)
            jj = kk <= Ntst && abs(τs[kk] - f1) < 1e-10 ? max(1, kk - 1) : searchsortedlast(τs, f1)
            clamp(jj, 1, Ntst)
        end
        kept = g <= jlast ? collect(g:jlast) : vcat(collect(g:Ntst), collect(1:jlast))
        isempty(kept) && error("Empty homoclinic time window [t0, t1]. Check the values of `t0`, `t1`.")
    end
    K = length(kept)
    @assert 1 <= K <= Ntst

    a0 = τs[g]      # arc start (boundary), physical fraction
    widths = [τs[idx + 1] - τs[idx] for idx in kept]
    L = sum(widths) # arc length, fraction of the period
    @assert L > 0

    # ━━━━━━ keep the stored points lying on the arc, in arc order ━━━━━━
    endlift = a0 + L
    idxs = Int[]
    lifts = Float64[]
    for i in eachindex(ft)
        li = ft[i] < a0 - 1e-12 ? ft[i] + 1 : ft[i]   # unwrapped fraction, starting at a0
        if a0 - 1e-12 <= li <= endlift + 1e-12
            push!(idxs, i)
            push!(lifts, li)
        end
    end
    perm = sortperm(lifts)
    idxs = idxs[perm]
    lifts = lifts[perm]
    # drop the duplicated node at the seam (fraction 0 ≡ fraction 1)
    uidx = Int[]
    ulifts = Float64[]
    for k in eachindex(idxs)
        if isempty(uidx) || lifts[k] - ulifts[end] > 1e-10
            push!(uidx, idxs[k])
            push!(ulifts, lifts[k])
        end
    end
    @assert length(uidx) == 1 + m * K "internal error: kept $(length(uidx)) nodes for $K collocation intervals"
    return (; kept, uidx, ulifts, widths, a0, L)
end

"""
$(TYPEDSIGNATURES)

Generate a [`HomoclinicHyperbolicProblemPBC`](@ref) functional together with an initial guess from a periodic orbit computed with the collocation method `coll`.

The guess is built by **extracting the time slices of the existing periodic-orbit solution**: a contiguous block of whole collocation intervals is selected and the corresponding stored points are reused as is — no polynomial re-interpolation is performed. A possibly *adapted* mesh (and therefore the resolution of peaked layers) is thus preserved.

!!! tip "Adapted mesh"
    In case of an adapted mesh, you can pass the `BK.POSavedSolutionAndState` solution directly in place of `x`, see the dedicated method below. When using the `AbstractArray` method, the mesh of `coll` must be the mesh on which `x` was computed; this can be forced with the keyword `mesh`.

## Arguments
- `coll::Collocation`: collocation discretization used to compute the periodic orbit
- `x::AbstractArray`: periodic orbit solution, as stored on a branch (e.g. `br.sol[end].x`)
- `pars`: parameters at which the periodic orbit was computed
- `lensHom::BK.AllOpticTypes`: parameter axis (lens) used for the continuation of the homoclinic orbit

## Keyword arguments
- `t0 = 0`, `t1 = 0`: absolute times in the periodic orbit delimiting the homoclinic interval. If both are `0`, the interval is obtained by removing the single collocation interval containing the point of the orbit closest to the saddle. Otherwise the whole collocation intervals covering the window `[t0, t1]` are kept.
- `ϵ0 = 1e-5`, `ϵ1 = 1e-5`: kept for compatibility; the distances `ϵ0hom`, `ϵ1hom` are measured from the first/last kept points to the saddle point.
- `maxT = Inf`: upper bound on the return time `T` of the homoclinic orbit
- `freeparams = ((@optic _.ϵ0), (@optic _.T))`: free parameters used to define the homoclinic orbit in parameter space
- `mesh = nothing`: fraction mesh (interval boundaries, length `Ntst + 1`) on which the solution `x` is defined. Automatically provided when `x` is a `BK.POSavedSolutionAndState`.
- `verbose = false`: print some debugging information

The extra `kwargs` are passed to the constructor of `::HomoclinicHyperbolicProblemPBC`.

## Output
- returns the tuple `(𝐇𝐨𝐦, xhom, pars, xhom)` where `𝐇𝐨𝐦::HomoclinicHyperbolicProblemPBC` and `xhom` is the initial guess. In the tutorials, only the first two entries are used.
"""
function generate_hom_problem(coll::Collocation,
                              x::AbstractArray,
                              pars,
                              lensHom::BK.AllOpticTypes;
                              verbose = false,
                              ϵ0 = 1e-5, ϵ1 = 1e-5,
                              t0 = 0, t1 = 0,
                              maxT = Inf,
                              freeparams = ((@optic _.ϵ0), (@optic _.T)),
                              mesh = nothing,
                              kw...)
    println("="^40)
    @assert coll.N > 0

    # working copy whose mesh is consistent with the solution `x`
    _coll = deepcopy(coll)
    if mesh !== nothing
        BK.update_mesh!(_coll, mesh)
    end

    T = getperiod(_coll, x)
    xc = BK.get_time_slices(_coll, x)          # n x (1 + m*Ntst)
    ft = BK.get_times(_coll)                   # node fractions, in [0, 1]
    τs = BK.getmesh(_coll)                     # interval boundaries
    n, m, Ntst = size(_coll)
    P = size(xc, 2)
    @assert P == 1 + m * Ntst
    # ━━━━━━ saddle point located on the existing nodes only ━━━━━━
    res = [norm(BK.residual(_coll.prob_vf, view(xc, :, i), pars)) for i in 1:P]
    isaddle = argmin(res)
    xs0 = xc[:, isaddle]
    # refine: exact equilibrium used to seed the saddle unknown and the projectors
    solS = BK._newton(_coll.prob_vf, xs0, pars, NewtonPar(;verbose))
    xsaddle = BK.converged(solS) ? solS.u : xs0
    # ━━━━━━ select the arc: contiguous whole intervals, extracted points only ━━━━━━
    arc = _select_homoclinic_arc(τs, ft, m, Ntst, isaddle, t0, t1, T)
    K = length(arc.kept)
    Nkeep = length(arc.uidx)
    a0 = arc.a0
    L = arc.L
    xflow = vec(xc[:, arc.uidx])              # extracted points, in arc order
    λ = (arc.ulifts .- a0) ./ L               # arc fractions of the kept nodes

    # reduced collocation discretization whose mesh follows the arc boundaries
    bvp = deepcopy(_coll)
    @reset bvp.update_section_every_step = 0
    bvp = BK._set_params_in_po(bvp, pars)
    if K != Ntst
        bvp = BK.set_collocation_size(bvp, K, m)
    end
    arc_mesh_fractions = zeros(eltype(_coll), K + 1)
    acc = zero(eltype(_coll))
    for (c, w) in enumerate(arc.widths)
        acc += w
        arc_mesh_fractions[c + 1] = acc / L
    end
    BK.update_mesh!(bvp, arc_mesh_fractions)

    # sanity check: the node times of `bvp` match the arc fractions of the extracted points
    @assert maximum(abs, λ .- BK.get_times(bvp)) < 1e-6 "node time mismatch in the homoclinic guess"

    Thom = min(L * T, maxT)
    x0 = xflow[1:n]
    x1 = xflow[end-n+1:end]

    BK.updatesection!(bvp, vcat(xflow, Thom), BK.getparams(bvp))

    # distances measured from the extracted endpoints
    ϵ0hom = norm(x0 - xsaddle)
    ϵ1hom = norm(x1 - xsaddle)

    # define problem for Homoclinic functional
    J = BK.jacobian(_coll.prob_vf, xsaddle, pars)
    𝐇𝐨𝐦 = HomoclinicHyperbolicProblemPBC(bvp,
                                          lensHom,
                                          length(xsaddle),
                                          copy(J);
                                          ϵ0 = ϵ0hom,
                                          ϵ1 = ϵ1hom,
                                          T = Thom,
                                          freeparams,
                                          kw...)

    @assert BK.getparams(𝐇𝐨𝐦) == pars "Errors with setting the parameters. Please an issue on the website of BifurcationKit."

    if verbose
        println("┌─ isaddle      = $isaddle")
        println("├─ kept         = $K intervals / $Nkeep nodes")
        println("├─ t0 (frac)    = $a0")
        println("├─ t1 (frac)    = $(mod(a0 + L, 1))")
        println("├─ T            = $Thom")
        println("└─ ϵ0hom, ϵ1hom = $((ϵ0hom, ϵ1hom))")
    end

    ns = 𝐇𝐨𝐦.nStable
    nu = 𝐇𝐨𝐦.nUnstable
    p1 = BK._get(pars, lensHom)

    xhom = ArrayPartition(xflow,
        xsaddle,
        zeros(eltype(xsaddle), n - ns, ns),
        zeros(eltype(xsaddle), n - nu, nu),
        [p1, map(x -> BK._get(𝐇𝐨𝐦, x), freeparams)...]
        )

    return 𝐇𝐨𝐦, xhom, pars, xhom
end

"""
$(TYPEDSIGNATURES)

Same as [`generate_hom_problem`](@ref) but `x` is a `BK.POSavedSolutionAndState`, as returned on the branch when mesh adaptation is used. The mesh `x._mesh` on which the solution `x.sol` is defined is restored in a working copy of `coll` and the section is updated with the phase `x.ϕ` before generating the homoclinic problem. The keyword arguments are the same as for the `AbstractArray` method.
"""
function generate_hom_problem(coll::Collocation,
                              x::BK.POSavedSolutionAndState,
                              pars,
                              lensHom::BK.AllOpticTypes;
                              k...)
    coll2 = deepcopy(coll)
    # enforce that the mesh of `coll2` matches the mesh on which `x.sol` is defined
    BK.update_mesh!(coll2, x._mesh)
    BK.updatesection!(coll2, x.ϕ, nothing)
    generate_hom_problem(coll2, x.sol, pars, lensHom; mesh = x._mesh, k...)
end
