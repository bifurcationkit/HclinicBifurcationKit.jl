# HclinicBifurcationKit.jl

This Julia package aims at performing **bifurcation analysis** of Homoclinic / Heteroclinic orbits of Cauchy problems.

It builds upon [BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/) (version ≥ 0.8) to perform continuation and numerical bifurcation analysis.

## 📦 Installation

Assuming that you already have Julia correctly installed, it suffices to install `HclinicBifurcationKit` in the standard way:

```julia
] add HclinicBifurcationKit
```

The package can also be installed directly from its repository:

```julia
] add https://github.com/bifurcationkit/HclinicBifurcationKit.jl
```

## Capabilities
- compute Homoclinic to Hyperbolic Saddle Orbits (HomHS) using orthogonal collocation or standard shooting
- compute bifurcations of HomHS
- start HomHS from a direct simulation
- automatic branch switching to HomHS from a Bogdanov–Takens bifurcation point

## 📚 Citing this work
If you use this package for your work, we ask that you **cite** the following paper! Open source development strongly depends on this. It is referenced on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with the *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib).

## 🧑‍💻 Other software

There are several good softwares already available.

- For continuation in small dimension, most software are listed on [DSWeb](https://dsweb.siam.org). One can mention the widely used AUTO-07p and [MATCONT](https://sourceforge.net/projects/matcont/). All these are very reliable and some address high co-dimension bifurcations.

- For large scale problems, there is none.

In Julia, the present package seems to be the only one.
