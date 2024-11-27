getF(hom::HomoclinicHyperbolicProblemPBC{<: ShootingProblem},x,p) = BK.vf(hom.bvp.flow, x, p)

################################################################################
"""
Custom section which align with the part of the orbit with largest norm. It thus keep track of the index of the time at which this occurs.
"""
mutable struct SectionSSmax{Tn, Tc}  <: BK.AbstractSection
    "Normal to define hyperplane"
    normal::Tn

    "Representative point on hyperplane"
    center::Tc

    "index of max norm"
    ind::Int
end

SectionSSmax(n,c) = SectionSSmax(n,c,0)
(sect::SectionSSmax)(u, T) = BK.sectionShooting(u, T, sect.normal, sect.center)
(sect::SectionSSmax)(u::AbstractMatrix, T) = BK.sectionShooting(u[:,sect.ind], T, sect.normal, sect.center)

# we update the field of Section, useful during continuation procedure for updating the section
function BK.update!(sect::SectionSSmax, normal, center, ind)
    copyto!(sect.normal, normal)
    copyto!(sect.center, center)
    sect.ind = ind
    sect
end

# this function updates the section during the continuation run
function BK.updatesection!(sh::ShootingProblem{Tf, Tjac, Ts, Tsection }, x, par) where {Tf <: BK.AbstractFlow, Tjac <: BK.AbstractJacobianType, Ts, Tsection <: SectionSSmax}
    xt = BK.get_time_slices(sh, x)
    ind = argmax(norm(xt[:, i]) for i=1:size(xt, 2))
    @views BK.update!(sh.section, BK.vf(sh.flow, xt[:, ind], par), xt[:, ind], ind)
    sh.section.normal ./= norm(sh.section.normal)
    return true
end

################################################################################
function get_homoclinic_orbit(hom::HomoclinicHyperbolicProblemPBC{Tbvp}, x::ArrayPartition, par0; kode...) where {Tbvp <: ShootingProblem}
    sh = hom.bvp
    M = BK.get_mesh_size(sh)

    # get the updated parameter
    param = x.x[5][1]
    lens = hom.lens
    par = set(par0, lens, param)

    xflow = x.x[1]
    N = div(length(xflow), M)
    T, ϵ0, ϵ1 = _changeHomParameters(hom, x.x[5])

    if  M>=1
        # return evolve(hom.bvp.flow, Val(:Full), xflow, par, T/M; kode...)
        return get_periodic_orbit(hom.bvp, vcat(xflow, T), par; kode...)
    else
        xshc = reshape(xflow, N, M)
        sol = [evolve(sh.flow, Val(:Full), xshc[:, ii], par, T/M; kode...) for ii in 1:M]
        time = sol[1].t; u = sol[1][:,:]
        for ii in 2:M
            append!(time, sol[ii].t .+ time[end])
            u = hcat(u, sol[ii][:,:])
        end
        return SolPeriodicOrbit(t = time, u = u)
    end
end

function initBVPforPBC(bvp0::ShootingProblem, prob_vf, Hom; N, T, ϵ)
    M = BK.get_mesh_size(bvp0)
    dt = 2T/M
    xflow = reduce(vcat, [Hom.orbit(-T + n*dt, ϵ) for n = 0:M-1] )
    bvp = @set bvp0.update_section_every_step = 0

    # section = isnothing(bvp.section) ? SectionSSmax(copy(xflow[1:N]), copy(xflow[1:N])) : bvp.section
    section = SectionSSmax(copy(xflow[1:N]), copy(xflow[1:N]))
    bvp = setproperties(bvp, section = section)
    BK.updatesection!(bvp, vcat(xflow, 2T), BK.getparams(bvp))
    return xflow, bvp
end

@views function (hom::HomoclinicHyperbolicProblemPBC{Tbvp, Nf})(x::ArrayPartition, par0) where {Tbvp <: ShootingProblem, Nf}
    @unpack N = hom
    sh = hom.bvp
    ns = hom.nStable
    nu = hom.nUnstable

    xflow = x.x[1]        # point x0 near the unstable manifold
    xsaddle = x.x[2]    # saddle point
    Ys = x.x[3]            #   stable part for CIS algo
    Yu = x.x[4]            # unstable part for CIS algo
    # get homoclinic parameters
    T, ϵ0, ϵ1 = _changeHomParameters(hom, x.x[5])
    @assert ϵ0 > 0 && ϵ1 > 0

    @assert size(Ys) == (N - ns, ns)
    @assert size(Yu) == (N - nu, nu)

    # get the updated parameter
    param = x.x[5][1]
    lens = hom.lens
    par = set(par0, lens, param)

    M = BK.get_mesh_size(sh)
    if M == 1
        x0 = xflow
        x1 = BK.evolve(sh.flow, x0, par, T).u # x1 = ϕ(x0, T)
    else
        N = div(length(xflow), M)
        @assert N*M == length(xflow)
        xshc = reshape(xflow, N, M)
        outshc = similar(xshc, N, M-1)
        if ~BK.isparallel(sh)
            for ii=1:M-1
                outshc[:, ii] .= BK.evolve(sh.flow, xshc[:, ii], par, T/M).u .- xshc[:, ii+1]
            end
            x1 = BK.evolve(sh.flow, xshc[:, M], par, T/M).u
        else
            solOde = BK.evolve(sh.flow, xshc, par, sh.ds .* T)
            for ii=1:M-1
                outshc[:, ii] .= solOde[ii][2] .- xshc[:, ii+1]
            end
            x1 = solOde[M][2]
        end
        x0 = xshc[:, 1]
    end

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
        if M==1
            outnrm[3] = sh.section(xflow, T)
        else
            # outnrm[3] = sh.section(BK.getTimeSlice(sh, xshc, 1), T)
            outnrm[3] = sh.section(xshc, 1)
        end
    end

    if M==1
        out = ArrayPartition(Fx, riU, riS, uP, sP, outnrm)
    else
        out = ArrayPartition(vec(outshc), Fx, riU, riS, uP, sP, outnrm)
    end
    return out
end

"""
$(SIGNATURES)

Generate a homoclinic to hyperbolic saddle problem from a periodic solution obtained with problem `pb`.

## Arguments
- `sh` a `ShootingProblem` which provide basic information, like the number of time slices `M`
- `x::AbstractArray` initial guess
- `pars` parameters
- `lensHom::BK.AllOpticTypes` parameter axis for continuation
- `ϵ0, ϵ1`: specify the distance to the saddle point of x₀, x₁
- `t0, t1`: specify the time corresponding to x₀, x₁. Overwrite the part with `ϵ0, ϵ1` if set.

## Optional arguments
You can pass the same arguments to the constructor of `::HomoclinicHyperbolicProblemPBC`.

## Output
- returns a `HomoclinicHyperbolicProblemPBC` and an initial guess.
"""
function generate_hom_problem(sh::ShootingProblem,
                            x::AbstractArray,
                            pars,
                            lensHom::BK.AllOpticTypes;
                            verbose = false,
                            time = LinRange(0, getperiod(sh, x), 100),
                            ϵ0 = 1e-5, ϵ1 = 1e-5,
                            t0 = 0, t1 = 0,
                            maxT = Inf,
                            freeparams = ((@optic _.ϵ0), (@optic _.T)),
                            kw...)
    verbose && println("="^40)
    @assert sh.M > 0
    T = getperiod(sh, x)
    M = BK.get_mesh_size(sh)

    # convert solution to homogenous mesh
    solpo = BK.get_po_solution(sh, x, pars)

    # find the saddle point as minimum of vector field norm
    xc = BK.get_time_slices(sh, x)
    ind_saddle = argmin(norm(BK.vf(sh.flow, solpo(t), pars)) for t in time)
    xsaddle = solpo(time[ind_saddle])
    tsaddle = time[ind_saddle]
    solsaddle = BK.newton(BifurcationProblem((x,p) -> BK.vf(sh.flow,x,p),xsaddle,pars), NewtonPar(verbose = true), norm = x->norm(x,Inf))
    if BK.converged(solsaddle)
        xsaddle .= solsaddle.u
    end

    if t1 == t0 == 0
        # find x0 and x1 on the unstable / stable subspace
        indUS = findfirst(norm(solpo(t) - xsaddle) > ϵ0 for t in time .+ tsaddle)
        t0 = mod(time[indUS]+tsaddle, T)
        x0 = solpo(t0)
        indS = findlast(norm(solpo(t) - xsaddle) > ϵ1 for t in time .+ t0)
        @assert ~isnothing(indS) "Increase ϵ0"
        t1 = time[indS]+t0
        x1 = solpo(t1)
    else
        x0 = solpo(t0)
        x1 = solpo(t1)
        indUS, indS = 0, 0
    end

    Thom = min(mod(t1-t0, T), maxT)

    # we put a uniform mesh in sh
    bvp = deepcopy(sh)
    bvp = BK.set_params_po(bvp, pars)
    @reset bvp.update_section_every_step = 0
    dt = Thom/M
    xflow = reduce(vcat, solpo(t0 + n * dt) for n=0:M-1)

    # put a specific section
    N = length(xsaddle)
    section = SectionSSmax(copy(xflow[1:N]), copy(xflow[1:N]))
    bvp = setproperties(bvp, section = section)
    BK.updatesection!(bvp, vcat(xflow, T), BK.getparams(bvp))

    # create Homoclinic parameters
    ϵ0hom = norm(x0 - xsaddle)
    ϵ1hom = norm(x1 - xsaddle)

    # define problem for Homoclinic functional
    J = ForwardDiff.jacobian(x -> BK.vf(sh.flow, x, pars), xsaddle)
    𝐇𝐨𝐦 = HomoclinicHyperbolicProblemPBC(bvp, lensHom, length(xsaddle), copy(J);  ϵ0 = ϵ0hom, ϵ1 = ϵ1hom, T = Thom, freeparams = freeparams, kw...)

    @assert BK.getparams(𝐇𝐨𝐦) == pars "Errors with setting the parameters. Please an issue on the website of BifurcationKit."

    if verbose
        println("┌─ tsaddle  = $tsaddle")
        println("├─ t0       = $t0")
        println("├─ t1       = $t1")
        println("├─ T        = $Thom")
        println("└─ is,i0,i1 = $((ind_saddle, indUS, indS))")
    end

    n = length(x0)
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
