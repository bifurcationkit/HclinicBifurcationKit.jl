# HclinicBifurcationKit.jl

This Julia package aims at performing **bifurcation analysis** of Homoclinic / Heteroclinic orbits of Cauchy problems.

It builds upon [BifurcationKit.jl]() with version > 0.2 to perform continuation and numerical bifurcation analysis.

## Installation

Assuming that you already have Julia correctly installed, it suffices to import  `HclinicBifurcationKit.jl` in the standard way:

`import Pkg; Pkg.add("https://github.com/bifurcationkit/HclinicBifurcationKit.jl")`

## Citing this work
If you use this package for your work, we ask that you **cite** the following paper!! Open source development strongly depends on this. It is referenced on HAL-Inria as follows:

```
@misc{veltz:hal-02902346,
  TITLE = {{BifurcationKit.jl}},
  AUTHOR = {Veltz, Romain},
  URL = {https://hal.archives-ouvertes.fr/hal-02902346},
  INSTITUTION = {{Inria Sophia-Antipolis}},
  YEAR = {2020},
  MONTH = Jul,
  KEYWORDS = {pseudo-arclength-continuation ; periodic-orbits ; floquet ; gpu ; bifurcation-diagram ; deflation ; newton-krylov},
  PDF = {https://hal.archives-ouvertes.fr/hal-02902346/file/354c9fb0d148262405609eed2cb7927818706f1f.tar.gz},
  HAL_ID = {hal-02902346},
  HAL_VERSION = {v1},
}
```

## Capabilities
- compute Homoclinic to Hyperbolic Saddle Orbits (HomHS) using Orthogonal collocation or Standard shooting
- compute bifurcation of HomHS
- start HomHS from a direct simulation
- automatic branch switching to HomHS from Bogdanov-Takes bifurcation point

## Other softwares

There are several good softwares already available.

- For continuation in small dimension, most softwares are listed on [DSWeb](https://ddebiftool.sourceforge.net). One can mention the widely used AUTO-07p and [MATCONT](https://sourceforge.net/projects/matcont/). All these are very reliable and some address high codimension bifurcations.

- For large scale problems, there is none.

In Julia, the present package seems to be the only one.