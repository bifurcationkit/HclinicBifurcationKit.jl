# Homoclinic to hyperbolic saddle

Consider the ODE problem written

$$\frac{du}{dt}=F(u(t),p)\tag{E}$$

where $p$ are parameters. 
A homoclinic solution $u^*$ to a hyperbolic saddle $u^s(p)$ satisfies $\lim\limits_{t\to\pm\infty}u^*(t) = u^s$ and $u^*(p)$ is a hyperbolic saddle of (E).

We provide 2 methods for computing such homoclinic orbits

2. one (Collocation) based on orthogonal collocation to discretize the above problem (E), with adaptive mesh 
3. one (Shooting) based on parallel standard shooting

## General method

The general method amounts to solving a boundary value problem which is simplified here for the exposition

$$\left\{\begin{aligned}
& \dot{u}(t)-2 T\cdot F(u(t), p)=0 \\
& F\left(u^s, p\right)=0 \\
& Q^{U^{\perp}, \mathrm{T}}\left(u(0)-u^s\right)=0, \\
& Q^{S^{\perp}, \mathrm{T}}\left(u(1)-u^s\right)=0 \\
& \left\|u(0)-u^s\right\|-\epsilon_0=0 \\
& \left\|u(1)-u^s\right\|-\epsilon_1=0 \\
\end{aligned}\right.$$

Basically, we truncate the homoclinic orbit on $[-T,T]$ and we impose that $u(-T),u(T)$ is close to $u^s$ and belong to the stable / unstable subspaces of $u^s$.
There are thus at most 3 free parameters and $T,\epsilon_0,\epsilon_1$ and the user can **either**

- chose one as a free parameter, for example $T$
- chose two as a free parameters, for example $T,\epsilon_1$



## Continuation

Please see the tutorials for examples. In a nutshell, you can compute homolinic orbits by setting up a [`HomoclinicHyperbolicProblemPBC`](@ref) or by branching from a Bogdanov-Takens point.

## Detection of codim 2 bifurcation points

You can detect the following codim 2 bifurcation points by using the option `detectCodim2Bifurcation` in the method `continuation`. Under the hood, the detection of these bifurcations is done by using Event detection as explained in [Event Handling](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/EventCallback/).

We refer to [^DeWitte] for a description of the bifurcations.


Type of bifurcation | Label
-------|-------------------
 Limit cycle |  LC 
 Homoclinic to Hyperbolic Saddle |  HHS 
 Homoclinic to Saddle-Node |  HSN 
 Neutral saddle |  NSS 
 Neutral saddle-focus |  NSF 
 Neutral Bi-Focus |  NFF 
 Shilnikov-Hopf |  SH 
 Double Real Stable leading eigenvalue |  DRS 
 Double Real Unstable leading eigenvalue |  DRU 
 Neutrally-Divergent saddle-focus (Stable) |  NDS 
 Neutrally-Divergent saddle-focus (Unstable) |  NDU 
 Three Leading eigenvalues (Stable) |  TLS 
 Three Leading eigenvalues (Unstable) |  TLU 
 Orbit-Flip with respect to the Stable manifold |  OFS 
 Orbit-Flip with respect to the Unstable manifold |  OFU 
 Non-Central Homoclinic to saddle-node |  NCH 

> Inclination-Flip with respect to the Stable / Unstable manifold is not yet detected.

## References

[^DeWitte]:> De Witte, Virginie, Willy Govaerts, Yuri A. Kuznetsov, and Mark Friedman. “Interactive Initialization and Continuation of Homoclinic and Heteroclinic Orbits in MATLAB.” ACM Transactions on Mathematical Software 38, no. 3 (April 2012): 1–34. https://doi.org/10.1145/2168773.2168776.