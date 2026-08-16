# HclinicBifurcationKit


| **Documentation** | **Build Status** | **Downloads** |
|:-----------------:|:----------------:|:-------------:|
| [![docs-dev][docs-dev-img]][docs-dev-url] |  [![Build Status](https://github.com/bifurcationkit/HclinicBifurcationKit.jl/workflows/CI/badge.svg)](https://github.com/bifurcationkit/HclinicBifurcationKit.jl/actions?query=workflow%3ACI) [![codecov](https://codecov.io/gh/bifurcationkit/HclinicBifurcationKit.jl/branch/main/graph/badge.svg?token=219HJEG8GM)](https://codecov.io/gh/bifurcationkit/HclinicBifurcationKit.jl) |  |

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://bifurcationkit.github.io/HclinicBifurcationKit.jl/dev

`HclinicBifurcationKit.jl` is a component package in the `BifurcationKit` ecosystem. It provides the utilities to compute homoclinic (heteroclinic) orbits of ODEs and to perform their numerical bifurcation analysis.

While completely independent and usable on its own, users interested in this functionality are encouraged to also check out [BifurcationKit.jl](https://github.com/bifurcationkit/BifurcationKit.jl).

## Features

- computation of Homoclinic to Hyperbolic Saddle Orbits (HomHS) using orthogonal collocation or (multiple) shooting
- bifurcation analysis of HomHS, including detection of codimension-two bifurcations
- automatic branch switching to HomHS from a Bogdanov–Takens point

## 📦 Installation

Assuming that you already have Julia correctly installed, it suffices to import
`HclinicBifurcationKit.jl` in the standard way:

```julia
] add HclinicBifurcationKit
```

## 📚 Support and citation
If you use `BifurcationKit.jl` in your work, we ask that you cite the following paper on [HAL-Inria](https://hal.archives-ouvertes.fr/hal-02902346) with *bibtex* entry [CITATION.bib](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/CITATION.bib). Open source development as part of academic research strongly depends on this. Please also consider starring this repository if you like our work, this will help us to secure funding in the future.
