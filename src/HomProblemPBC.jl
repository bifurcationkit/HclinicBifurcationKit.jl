"""
Computation of homoclinic orbit to an hyperbolic saddle based on the projection boundary condition (PBC) method.
$(SIGNATURES)

$(TYPEDFIELDS)
"""
mutable struct HomoclinicHyperbolicProblemPBC{Tbvp, Nfree, Tlens, Ty, Tlensfree, Tq, Tt} <: BK.AbstractBoundaryValueProblem
    "Sructure encoding the boundary value problem. For example, you can pass a `PeriodicOrbitTrapProblem`, a `PeriodicOrbitOCollProblem` or an `AbstractShootingProblem`"
    bvp::Tbvp

    "Two lenses which are used to define 2 free parameters."
    lens::Tlens

    "Return time T"
    T::Ty

    "Precision of how far the section is from the homoclinic point"
    ϵ0::Ty

    "Precision of how far the section is from the homoclinic point"
    ϵ1::Ty

    "Free parameters"
    freelens::Tlensfree

    "Orthonormal Projector on the unstable subspace orthogonal"
    Qu0::Tq # size = n x ns

    "Orthonormal Projector on the stable subspace orthogonal"
    Qs0::Tq # size = n x nu

    "Dimension of phase space"
    N::Int64

    "updates the section every `update_section_every_step` step during continuation"
    updateEveryStep::Int

    "How the jacobian of the problem is computed"
    jacobian::Symbol

    test::Tt
    testOrbitFlip::Bool
    testInclinationFlip::Bool
    nUnstable::Int64
    nStable::Int64
end
@inline BK.getparams(pb::HomoclinicHyperbolicProblemPBC) = BK.getparams(pb.bvp)
@inline BK.getlens(pb::HomoclinicHyperbolicProblemPBC) = pb.lens

function HomoclinicHyperbolicProblemPBC(bvp::Tbvp,
                lens::Tlens,
                N::Int,
                J::AbstractMatrix{Ty}; ϵ0 = 0.01,
                ϵ1 = 0.01,
                T = 100.,
                freeparams = ((@lens _.ϵ0), (@lens _.T)),
                update_every_step = 2,
                testOrbitFlip = false,
                testInclinationFlip = false,
                jacobian::Symbol = :autodiffDense
                ) where {Ty, Tbvp, Tlens}
    @assert ~(bvp isa PeriodicOrbitTrapProblem) "This type of BVP is not handled yet."
    T = convert(Ty, T)
    ϵ0 = convert(Ty, ϵ0)
    ϵ1 = convert(Ty, ϵ1)
    @unpack Qs0, Qu0, nStable, nUnstable = get_S_U_stableSpaces(J)
    Nfree = length(freeparams)
    @assert Nfree < 3 "At most 2 homoclinic parameters"
    test0 = (NNS=one(T), NSF=one(T), NFF=one(T), DRS=one(T), DRU=one(T), NDS=one(T), NDU=one(T), TLS=one(T), TLU=one(T), NCH=one(T), SH=one(T), BT=one(T),OFU=one(T), OFS=one(T), IFU=one(T), IFS=one(T))
    HomoclinicHyperbolicProblemPBC{Tbvp, Nfree, Tlens, Ty, typeof(freeparams), typeof(Qs0), typeof(test0)}(bvp, lens, T, ϵ0, ϵ1, freeparams, Qu0, Qs0, N, update_every_step, jacobian, test0,testOrbitFlip, testInclinationFlip, nUnstable, nStable )
end

function Base.show(io::IO, hom::HomoclinicHyperbolicProblemPBC)
    println(io, "┌─ Homoclinic problem (PBC)")
    println(io, "├─ ϵ0          : ", hom.ϵ0)
    println(io, "├─ ϵ1          : ", hom.ϵ1)
    println(io, "├─ T           : ", hom.T)
    println(io, "├─ free params : ", map(BK.get_lens_symbol,hom.freelens))
    println(io, "├─ update prob : ", hom.updateEveryStep)
    println(io, "└─ lens        : ", BK.get_lens_symbol(BK.getlens(hom)))
    printstyled(io, "Boundary value problem:\n"; color = :blue, bold = true)
    show(io, hom.bvp)
end

"""
Compute projectors on stable / unstable linear manifolds for the jacobian `J`.
"""
function get_S_U_stableSpaces(J)
    F = schur(J)
    @assert abs(prod(F.values))>0 "This is not a saddle point"
    Qu0, nUnstable = _getSUSpace(F, >)
    Qs0, nStable   = _getSUSpace(F, <)
    return (; Qs0, Qu0, nStable, nUnstable)
end

function _getSUSpace(F, op)
    sel = op.(real.(F.values), 0)
    n = sum(sel)
    Fp = ordschur(F, sel)
    return Fp.vectors, n
end

# we follow De Witte, Virginie, Willy Govaerts, Yuri A. Kuznetsov, and Mark Friedman. “Interactive Initialization and Continuation of Homoclinic and Heteroclinic Orbits in MATLAB.” ACM Transactions on Mathematical Software 38, no. 3 (April 2012): 1–34. https://doi.org/10.1145/2168773.2168776.

@views function ricattiBlocks(Q, J, n)
    #TODO  Faire  inplace
    T = Q' * J * Q
    T11  = T[1:n, 1:n]
    T12  = T[1:n, n+1:end]
    T21  = T[n+1:end, 1:n]
    T22  = T[n+1:end, n+1:end]
    return (; T11, T12, T21, T22)
end

function ricattiEq(Ts, Y)
    T11, T12, T21, T22 = Ts
    return T22 * Y - Y * T11 + T21 - Y * T12 * Y
end

function _changeHomParameters(hom::HomoclinicHyperbolicProblemPBC, xpar)
    @unpack T, ϵ0, ϵ1 = hom
    lensp = hom.freelens
    for (id, l) in enumerate(lensp)
        T, ϵ0, ϵ1 = set((;T, ϵ0, ϵ1), l, xpar[1+id])
    end
    return (T  = T, ϵ0 = abs(ϵ0), ϵ1 = abs(ϵ1))
end

function generate_hom_solution(pb::BK.AbstractBoundaryValueProblem, orbit, T)
    M = BK.get_mesh_size(pb)
    orbitguess_a = [orbit(-T + t * (2T)) for t in LinRange(0, 1, M + 1)[1:M]]
    # append period at the end of the initial guess
    orbitguess_v = reduce(vcat, orbitguess_a)
    if pb isa BK.AbstractPoincareShootingProblem
        return vec(orbitguess_v)
    end
    orbitguess_v
end
####################################################################################################

getVectorField(bvp::ShootingProblem,x,p) = BK.vf(bvp.flow,x,p)
getVectorField(bvp::PeriodicOrbitOCollProblem,x,p) = BK.residual(bvp.prob_vf,x,p)

function get_tests_for_HHS(𝐇𝐨𝐦, J, z, xₛ, x₀, x₁, T::Ty, pars; tol = 1e-5) where Ty
    F = eigen(J)
    ind = sortperm(F.values, by = real)
    values = F.values[ind]
    ζs = F.vectors[:, ind]
    ind0 = findfirst( real(v) > 0 for v in values )
    @assert ind0>1 "We need a stable projector! The point is not a saddle"
    inf = typemax(Inf)
    λ₁ = values[ind0]; ζ₁ᵘ = ζs[:, ind0]
    λ₂ = ind0 == length(values) ? inf : values[ind0+1]
    λ₃ = ind0 >= length(values)-1 ? inf : values[ind0+2]
    μ₁ = values[ind0-1]; ; ζ₁ˢ = ζs[:, ind0-1]
    μ₂ = ind0 <= 2 ? -inf : values[ind0-2]
    μ₃ = ind0 <= 3 ? -inf : values[ind0-3]

    # adjoint case
    F★ = eigen(J')
    ind★ = sortperm(F★.values, by = real)
    values★ = F★.values[ind★]
    ζ★s = F★.vectors[:, ind★]
    # display(hcat(values, values★))
    ζ★₁ᵘ = ζ★s[:, ind0]
    ζ★₁ˢ = ζ★s[:, ind0-1]

    @assert norm(values★ - values, Inf) < 1e-10

    normalize!(ζ₁ᵘ);
    normalize!(ζ₁ˢ);
    ζ★₁ᵘ ./= dot(ζ★₁ᵘ, ζ₁ᵘ)
    ζ★₁ˢ ./= dot(ζ★₁ˢ, ζ₁ˢ)

    @assert values★[ind0] ≈ λ₁
    @assert values★[ind0-1] ≈ μ₁

    NNS = real(μ₁) + real(λ₁)
    if abs(imag(μ₁)) < tol && abs(imag(λ₁)) < tol
        NNS = real(μ₁) + real(λ₁)
        NSF = -inf
        NFF = -inf
    elseif abs(imag(μ₁)) > tol && abs(imag(λ₁)) < tol
        NNS = -inf
        NSF = real(μ₁) + real(λ₁)
        NFF = -inf
    elseif abs(imag(μ₁)) < tol && abs(imag(λ₁)) > tol
        NNS = -inf
        NSF = real(μ₁) + real(λ₁)
        NFF = -inf
    else
        NNS = -inf
        NSF = -inf
        NFF = real(μ₁) + real(λ₁)
    end
    DRS = (abs(imag(μ₁)) < tol) ? (real(μ₁) - real(μ₂))^2 : -(imag(μ₁) - imag(μ₂))^2
    DRU = (abs(imag(λ₁)) < tol) ? (real(λ₁) - real(λ₂))^2 : -(imag(λ₁) - imag(λ₂))^2
    NDS = real(μ₁) + real(μ₂) + real(λ₁)
    NDU = real(μ₁) + real(λ₂) + real(λ₁)
    TLS = real(μ₁) - real(μ₃)
    TLU = real(λ₁) + real(λ₃)
    NCH = real(μ₁)
    SH  = real(λ₁)
    BT  = real(μ₁) * real(λ₁)

    # orbit flip
    OFS = 1
    OFU = 1

    if 𝐇𝐨𝐦.testOrbitFlip
        _ii = maximum(real(v) for v in ζ★₁ˢ)

        OFS = dot(real(ζ★₁ˢ), x₁ - xₛ) #* exp(-real(μ₁) * T)
        if _ii < 0
            OFS *= -1
        end

        _ii = argmax(abs(v) for v in ζ★₁ᵘ)
        # ζ★₁ᵘ .*= sign(real(ζ★₁ᵘ[_ii]))
        OFU = dot(real(ζ★₁ᵘ), x₀ - xₛ) #* exp(real(λ₁) * T)
        if abs(imag(λ₁)) > tol
            OFU *= dot(imag(ζ★₁ᵘ), x₀ - xₛ) #* exp(real(λ₁) * T)
        end
    end

    # inclination flip
    IFU = 1
    IFS = 1

    𝐇𝐨𝐦.test = (;NNS, NSF, NFF, DRS, DRU, NDS, NDU, TLS, TLU, NCH, SH, BT, OFU, OFS, IFU, IFS)
    # for (n,v) in pairs(𝐇𝐨𝐦.test)
    #     println(n, " - ", v)
    # end
    return 𝐇𝐨𝐦.test
end

"""
$(SIGNATURES)

This is the continuation method for computing a homoclinic solution to a hyperbolic saddle. The parameter lens is the one from `prob_vf::BifurcationProblem`.

# Arguments

Similar to [`continuation`](@ref) except that the problem is a [`HomoclinicHyperbolicProblemPBC`](@ref).
"""
function BK.continuation(𝐇𝐨𝐦::HomoclinicHyperbolicProblemPBC,
            homguess,
            lens::Lens,
            alg::BK.AbstractContinuationAlgorithm,
            _contParams::ContinuationPar;
            plot_solution = BK.plot_default,
            kwargs...
            )
    function updateHom(z, tau, step, contResult; kUP...)
        # if this is called from bisection, do not update the problem
        bisection = get(kUP, :bisection, false)
        if bisection
            return true
        end
        verbose = get(kwargs, :verbosity, 0) >= 1
        # user-passed finalizer
        finaliseUser = get(kwargs, :finaliseSolution, nothing)

        # we first check that the continuation step was successful
        # if not, we do not update the problem with bad information!
        success = get(kUP, :state, nothing).converged
        if (~BK.mod_counter(step, 𝐇𝐨𝐦.updateEveryStep) || success == false)
            return isnothing(finaliseUser) ? true : finaliseUser(z, tau, step, contResult; prob = 𝐇𝐨𝐦, kUP...)
        end

        T, ϵ0, ϵ1 = _changeHomParameters(𝐇𝐨𝐦, z.u.x[5])
        @assert (ϵ0 > 0) && (ϵ1 > 0)
        newpar = set(BK.getparams(𝐇𝐨𝐦), BK.getlens(𝐇𝐨𝐦), z.u.x[end][1])
        newpar = set(newpar, lens, z.p)
        verbose && (@info "Update 𝐇𝐨𝐦" T, ϵ0, ϵ1)

        # compute the jacobian
        J = ForwardDiff.jacobian(z -> getVectorField(𝐇𝐨𝐦.bvp, z, newpar ), z.u.x[2])
        @unpack Qs0, Qu0, nStable, nUnstable = get_S_U_stableSpaces(J)

        # this is a Hack for Orthogonal collocation
        if success && (𝐇𝐨𝐦.bvp isa PeriodicOrbitOCollProblem) && 𝐇𝐨𝐦.bvp.meshadapt
            verbose && (@info "update mesh!")
            oldsol = BK._copy(z)
            oldmesh = BK.get_times(𝐇𝐨𝐦.bvp) .* getperiod(𝐇𝐨𝐦.bvp, z.u, nothing)
            oldu = vcat(z.u.x[1], T)
            adapt = BK.compute_error!(𝐇𝐨𝐦.bvp, oldu;
                    verbosity = 𝐇𝐨𝐦.bvp.verbose_mesh_adapt,
                    par = newpar,
                    K = 𝐇𝐨𝐦.bvp.K)
            z.u.x[1] .= oldu[1:end-1]
            if ~adapt.success
                return false
            end
            @info norm(oldsol.u - z.u, Inf)
        end

        if success && BK.mod_counter(step, 𝐇𝐨𝐦.updateEveryStep) == 1
            𝐇𝐨𝐦.Qs0 .= Qs0
            𝐇𝐨𝐦.Qu0 .= Qu0
            𝐇𝐨𝐦.nUnstable = nUnstable
            𝐇𝐨𝐦.nStable = nStable
            𝐇𝐨𝐦.T = T
            𝐇𝐨𝐦.ϵ0 = ϵ0
            𝐇𝐨𝐦.ϵ1 = ϵ1

            # update section
            BK.updatesection!(𝐇𝐨𝐦.bvp, vcat(z.u.x[1], T), newpar)
        end

        # call the user-passed finalizer
        resFinal = isnothing(finaliseUser) ? true : finaliseUser(z, tau, step, contResult; prob = 𝐇𝐨𝐦, kUP...)
    end

    function testHom(iter, state)
        z = getx(state)

        T, ϵ0, ϵ1 = _changeHomParameters(𝐇𝐨𝐦, z.x[5])
        newpar = set(BK.getparams(𝐇𝐨𝐦), BK.getlens(𝐇𝐨𝐦), z.x[end][1])
        newpar = set(newpar, lens, getp(state))

        # compute the jacobian
        J = ForwardDiff.jacobian(z -> getVectorField(𝐇𝐨𝐦.bvp, z, newpar ), z.x[2])

        u = vcat(z.x[1], T)
        uc = BK.get_time_slices(𝐇𝐨𝐦.bvp, u)
        x0 = @view uc[:, 1]
        x1 = @view uc[:, end]

        get_tests_for_HHS(𝐇𝐨𝐦, J, z, z.x[2], x0, x1, T, newpar) |> values
    end

    probhom_bk = BifurcationProblem(𝐇𝐨𝐦, homguess, BK.getparams(𝐇𝐨𝐦), lens;
        J = (x,p) -> ForwardDiff.jacobian(z -> 𝐇𝐨𝐦(z,p), x),
        record_from_solution = (x, p) -> begin
            if length(𝐇𝐨𝐦.freelens) == 1
                lensS = map(BK.get_lens_symbol, (BK.getlens(𝐇𝐨𝐦), lens, @lens _.FreeP1))
            else
                lensS = map(BK.get_lens_symbol, (BK.getlens(𝐇𝐨𝐦), lens, (@lens _.FreeP1), @lens _.FreeP2))
            end
            record = (;zip(lensS, (x.x[end][1], p, x.x[end][2:end]...))...)
            if _contParams.detect_event > 0
                record = merge(record, 𝐇𝐨𝐦.test)
            end
            return record
        end,
        plot_solution = modify_hom_plot(𝐇𝐨𝐦, lens, (kwargs..., plot_solution = plot_solution)),
        )

    event = ContinuousEvent(16, testHom, false, ("NNS", "NSF", "NFF", "DRS", "DRU", "NDS", "NDU", "TLS", "TLU", "NCH", "SH", "BT", "OFU", "OFS", "IFU", "IFS"), 0)

    return BK.continuation(probhom_bk, alg, _contParams;
        finalise_solution = 𝐇𝐨𝐦.updateEveryStep == 0 ? get(kwargs, :finalise_solution, BK.finalise_default) : updateHom,
        kind = HomoclinicHyperbolicSaddleCont(),
        event = event,
        kwargs...)
end
####################################################################################################
"""
$(SIGNATURES)

Perform automatic branch switching to homoclinic curve from a Bogdanov-Takens bifurcation point. It uses the homoclinic orbit predictor from the Bogdanov-Takens normal form.

# Arguments
- `prob::BifurcationProblem` contains the vector field
- `bt::BK.BogdanovTakens` a Bogdanov-takens point. For example, you can get this from a call to `bt = get_normal_form(br, ind_bt)`
- `bvp::BK.AbstractBoundaryValueProblem`, for example `PeriodicOrbitOCollProblem(50, 4)`
- `alg` continuation algorithm
- `_contParams::ContinuationPar`

## Optional arguments
- `ϵ0 = 1e-5` distance of the homolinic orbit from the saddle point
- `amplitude = 1e-3` amplitude of the homoclinic orbit
- `maxT = Inf` limit on the "period" of the homoclinic cycle.

You can also pass the same arguments to the constructor of `::HomoclinicHyperbolicProblemPBC` and those to `continuation` from BifurcationKit.
- `kwargs` arguments passed to `continuation`
"""
function BK.continuation(prob_vf,
            bt::BK.BogdanovTakens,
            bvp::BK.AbstractBoundaryValueProblem,
            alg::BK.AbstractContinuationAlgorithm,
            _contParams::ContinuationPar ;
            ϵ0 = 1e-5, amplitude = 1e-3,
            freeparams = ((@lens _.ϵ0), (@lens _.T)),
            maxT = Inf,
            update_every_step = 1,
            test_orbit_flip = false,
            test_inclination_flip = false,
            kwargs...
            )
    printstyled(color=:magenta, "\n\n────────────────────────────\n┌─ Debut Hom_BT init\n")

    Hom = predictor(bt, Val(:HomoclinicCurve), 0.)
    ϵ = sqrt(amplitude * abs(bt.nf.a) / 6)
    printstyled(color = :magenta,"├─ ϵ = ", ϵ, "\n")

    # get the homoclinic parameters
    p1, p2 = Hom.α(ϵ)
    lens1, lens2 = bt.lens
    pars = set(bt.params, lens1, p1)
    pars = set(pars, lens2, p2)

    # estimate half return time T by solving sech(εT)= √(κ/amp)
    # sech(x) ~ 2/exp(x)
    κ = ϵ0
    ζ = sqrt(κ / amplitude)
    printstyled(color = :magenta,"├───  ζ = ", ζ, "\n")
    T0 = 2/exp(ζ)
    f(t) = sech(ϵ * t) - ζ
    pb = BifurcationProblem((t,p) -> [f(t[1])], [T0/ϵ], nothing)
    solT = newton(pb, NewtonPar(tol = 1e-8, verbose = false))
    @assert BK.converged(solT) "Newton iteration failed in determining the half return time T"
    T = min(solT.u[1], maxT)
    printstyled(color = :magenta,"├─ T = ", T, "\n")
    T < 40 && (@warn "Homoclinic cycle approximation is very crude! Decrease ϵ0")

    # get the saddle point
    xsaddle = Hom.orbit(100_000., ϵ)
    N = length(xsaddle)

    # get the points x(0) and x(1)
    x0 = Hom.orbit(-T, ϵ)
    x1 = Hom.orbit( T, ϵ)

    # create Homoclinic parameters
    ϵ0hom = norm(x0 - xsaddle)
    ϵ1hom = norm(x1 - xsaddle)
    Thom = 2T

    # define initial guess for newton
    prob_vf = re_make(prob_vf, params = pars)
    xflow, bvp = initBVPforPBC(bvp, prob_vf, Hom; N = N, T = T, ϵ = ϵ)
    bvp = BK.set_params_po(bvp, pars)

    # define problem for Homoclinic functional
    J = BK.jacobian(prob_vf, xsaddle, pars)
    𝐇𝐨𝐦 = HomoclinicHyperbolicProblemPBC(bvp, lens1, length(xsaddle), copy(J);  ϵ0 = ϵ0hom, ϵ1 = ϵ1hom, T = Thom, freeparams = freeparams, update_every_step = update_every_step, testOrbitFlip = test_orbit_flip,
    testInclinationFlip = test_inclination_flip)

    @assert BK.getparams(𝐇𝐨𝐦) == pars "Errors with setting the parameters. Please an issue on the website of BifurcationKit."

    printstyled("──> convergence to saddle point:\n", color = :magenta)
    solsaddle = BK.newton(BifurcationProblem((x,p) -> getF(𝐇𝐨𝐦,x,p), xsaddle, pars), NewtonPar(verbose = true), norm = x->norm(x,Inf))
    if BK.converged(solsaddle)
        xsaddle .= solsaddle.u
    end

    ns = 𝐇𝐨𝐦.nStable
    nu = 𝐇𝐨𝐦.nUnstable

    xhom = ArrayPartition(xflow,
        xsaddle,
        zeros(eltype(x0), N - ns, ns),
        zeros(eltype(x0), N - nu, nu),
        [p1, map(x -> get(𝐇𝐨𝐦,x), freeparams)...]
        )

    br = BK.continuation(𝐇𝐨𝐦, xhom, lens2, alg, _contParams; kwargs...)
    return Branch(br, bt)
end
