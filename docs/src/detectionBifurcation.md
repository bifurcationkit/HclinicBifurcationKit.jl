# Detection of bifurcation points

The bifurcations are detected during a call to `br = continuation(prob, alg, contParams::ContinuationPar;kwargs...)` by turning on the following flags:

- `contParams.detectBifurcation = 2` (for eigenvalues based bifurcations)
- `contParams.detectEvent = 2` (for other bifurcations like inclinations)
    
## Precise detection of bifurcation points using Bisection    

Note that the bifurcation points detected when `detectBifurcation = 2` can be rather *crude*  localization of the true bifurcation points. Indeed, we only signal that, in between two continuation steps *which can be large*, a (several) bifurcation has been detected. Hence, we only have a rough idea of where the bifurcation is located, unless your `dsmax` is very small... This can be improved as follows.

If you choose `detectBifurcation = 3`, a **bisection algorithm** is used to locate the bifurcation points more precisely. It means that we recursively track down the change in stability. Some options in [`ContinuationPar`](@ref) control this behavior:

- `nInversion`: number of sign inversions in the bisection algorithm
- `maxBisectionSteps` maximum number of bisection steps
- `tolBisectionEigenvalue` tolerance on real part of eigenvalue to detect bifurcation points in the bisection steps

!!! tip "Bisection mode"
    During the bisection, the eigensolvers are called like `eil(J, nev; bisection = true)` in order to be able to adapt the solver precision.


## List of detected bifurcation points

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