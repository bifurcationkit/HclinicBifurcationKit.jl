# HclinicBifurcationKit


| **Documentation** | **Build Status** | **Downloads** |
|:-----------------:|:----------------:|:-------------:|
| [![docs-dev][docs-dev-img]][docs-dev-url] |  [![Build Status](https://github.com/bifurcationkit/HclinicBifurcationKit.jl/workflows/CI/badge.svg)](https://github.com/bifurcationkit/HclinicBifurcationKit.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/bifurcationkit/HclinicBifurcationKit.jl/branch/main/graph/badge.svg?token=219HJEG8GM)](https://codecov.io/gh/bifurcationkit/HclinicBifurcationKit.jl) |  |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://bifurcationkit.github.io/HclinicBifurcationKit.jl/dev


`HclinicBifurcationKit.jl` is a component package in the `BifurcationKit` ecosystem. It holds the utilities concerning the computations of homoclinic (heteroclinic) orbits. While completely independent
and usable on its own, users interested in using this
functionality should check out [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl).

## Installation

Assuming that you already have Julia correctly installed, it suffices to import
`HclinicBifurcationKit.jl` in the standard way:

```julia
] add HclinicBifurcationKit
```

## Support and citation
If you use this package for your work, we ask that you cite the following paper. Open source development as part of academic research strongly depends on this. Please also consider starring this repository if you like our work, this will help us to secure funding in the future. It is referenced on HAL-Inria as follows:

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
