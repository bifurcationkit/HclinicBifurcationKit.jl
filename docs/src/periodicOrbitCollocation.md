# Homoclinic based on orthogonal collocation

We compute homoclinic orbits by discretizing the Cauchy problem on `Ntst` time intervals with orthogonal collocation, as implemented in the structure `BifurcationKit.Collocation`.

!!! warning "Large scale"
    The current implementation is not yet optimized for large scale problems. This will be improved in the future.

The general method is explained in the [periodic orbit collocation](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/periodicOrbitCollocation/) section of the BifurcationKit.jl documentation.

## General method

Please see [^DeWitte] for a thorough description of the method. It amounts to solving a boundary value problem.

$$\left\{\begin{aligned}
& \dot{u}(t)-2 T\cdot F(u(t), p)=0 \\
& F\left(u_0, p\right)=0 \\
& Q^{U^{\perp}, \mathrm{T}}\left(u(0)-u_0\right)=0, \\
& Q^{S^{\perp}, \mathrm{T}}\left(u(1)-u_0\right)=0 \\
& T_{22 U} Y_U-Y_U T_{11 U}+T_{21 U}-Y_U T_{12 U} Y_U=0, \\
& T_{22 S} Y_S-Y_S T_{11 S}+T_{21 S}-Y_S T_{12 S} Y_S=0 \\
& \left\|u(0)-u_0\right\|-\epsilon_0=0 \\
& \left\|u(1)-u_0\right\|-\epsilon_1=0 \\
& \int_0^1 \tilde{u}^*(t)[u(t)-\tilde{u}(t)] d t=0, \\
\end{aligned}\right.$$

## Mesh adaptation

The goal of this functionality is to adapt the mesh in order to minimize the discretization error. It is activated with `meshadapt = true` when constructing a `BifurcationKit.Collocation`; the number of allowed adaptation steps is controlled by the parameter `K` of the collocation.

When mesh adaptation is used, the solutions stored on a branch are `BifurcationKit.POSavedSolutionAndState` which record the mesh and the phase condition. [`generate_hom_problem`](@ref) accepts such a solution directly, so that the adapted mesh and phase are properly restored before building the homoclinic problem. The initial guess is built by **extracting the time slices of the periodic orbit** (no polynomial re-interpolation): a contiguous block of whole collocation intervals is selected, thereby preserving the resolution of peaked layers provided by the adaptive mesh. When no time window is given, the single collocation interval containing the point closest to the saddle point is removed; otherwise the whole intervals covering the window `[t0, t1]` are kept.

## Usage

A typical workflow is to

1. compute a branch of periodic orbits with `BifurcationKit` using a `Collocation` discretization,
2. build the homoclinic problem from the last point of the branch with [`generate_hom_problem`](@ref),
3. continue the homoclinic orbit with [`continuation`](@ref):

```julia
# coll is a Collocation discretization, brpo a branch of periodic orbits
𝐇𝐨𝐦, xhom, pars, _ = generate_hom_problem(coll, brpo.sol[end].x, BK.setparam(brpo, brpo.sol[end].p), BK.getlens(brpo))
br_hom = continuation(𝐇𝐨𝐦, xhom, lens, PALC(), ContinuationPar(); kwargs...)
```

See the [tutorials](@ref tutorials-page) for fully worked examples.

## Jacobian

The jacobian is computed with automatic differentiation *e.g.* `ForwardDiff.jl`

## References

[^DeWitte]:> De Witte, Virginie, Willy Govaerts, Yuri A. Kuznetsov, and Mark Friedman. “Interactive Initialization and Continuation of Homoclinic and Heteroclinic Orbits in MATLAB.” ACM Transactions on Mathematical Software 38, no. 3 (April 2012): 1–34. https://doi.org/10.1145/2168773.2168776.
