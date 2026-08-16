# Homoclinic to hyperbolic saddle

Consider the ODE problem written

$$\frac{du}{dt}=F(u(t),p)\tag{E}$$

where $p$ denotes the parameters.
A homoclinic solution $u^*$ to a hyperbolic saddle $u^s(p)$ satisfies $\lim\limits_{t\to\pm\infty}u^*(t) = u^s$ and $u^*(p)$ is a hyperbolic saddle of (E).

We provide 2 methods for computing such homoclinic orbits:

1. [Homoclinic based on orthogonal collocation](@ref): orthogonal collocation to discretize the above problem (E), with adaptive mesh,
2. [Homoclinic based on parallel multiple shooting](@ref): standard shooting based on the flow of (E).

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

Basically, we truncate the homoclinic orbit on $[-T,T]$ and we impose that $u(-T)$ and $u(T)$ are close to $u^s$ and belong to the stable / unstable subspaces of $u^s$.

The homoclinic solution is thus parametrized by the three scalars $T$, $\epsilon_0$ and $\epsilon_1$. Besides the continuation parameter, the user must select the free parameters of the problem (keyword `freeparams` in [`HomoclinicHyperbolicProblemPBC`](@ref)); at most two of $T, \epsilon_0, \epsilon_1$ can be free, e.g.

- one free parameter, for example $T$,
- two free parameters, for example $T,\epsilon_1$.

## Continuation

The homoclinic problem can be set up with [`generate_hom_problem`](@ref) (from a periodic orbit) or by branch switching from a Bogdanov–Takens point with [`continuation`](@ref); see the [tutorials](@ref tutorials-page) for examples.

## Detection of codimension 2 bifurcation points

Codimension-two bifurcation points along the homoclinic branch can be detected during the continuation. We refer to the page [Detection of bifurcation points](@ref) for the list of detected bifurcations and how to enable them.

## References

[^DeWitte]:> De Witte, Virginie, Willy Govaerts, Yuri A. Kuznetsov, and Mark Friedman. “Interactive Initialization and Continuation of Homoclinic and Heteroclinic Orbits in MATLAB.” ACM Transactions on Mathematical Software 38, no. 3 (April 2012): 1–34. https://doi.org/10.1145/2168773.2168776.
