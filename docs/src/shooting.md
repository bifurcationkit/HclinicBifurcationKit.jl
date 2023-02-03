# Homoclinic based on parallel multiple shooting

We aim at finding homoclinic orbits for the Cauchy problem 

$$\tag{1} \frac{d x}{d t}=f(x)$$ 

and we write $\phi^t(x_0)$ the associated flow (or semigroup of solutions).

!!! warning "Large scale"
    The current implementation is not yet optimised for large scale problems. This will be improved in the future.    

The general method is explained in [BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/periodicOrbitShooting/).

## General method

It amounts to solving a boundary value problem. See [^DeWitte] for description of the equations on the projectors.

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

## Jacobian

The jacobian is computed with automatic differentiation *e.g.* `ForwardDiff.jl`

## References

[^DeWitte]:> De Witte, Virginie, Willy Govaerts, Yuri A. Kuznetsov, and Mark Friedman. “Interactive Initialization and Continuation of Homoclinic and Heteroclinic Orbits in MATLAB.” ACM Transactions on Mathematical Software 38, no. 3 (April 2012): 1–34. https://doi.org/10.1145/2168773.2168776.
