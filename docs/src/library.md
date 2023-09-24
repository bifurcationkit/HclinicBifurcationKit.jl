# Library

```@contents
Pages = ["library.md"]
Depth = 3
```

## Parameters


## Problems

```@docs
HclinicBifurcationKit.HomoclinicHyperbolicProblemPBC
```

## Continuation

```@docs
BifurcationKit.continuation(𝐇𝐨𝐦::HclinicBifurcationKit.HomoclinicHyperbolicProblemPBC,homguess,lens::Lens,alg::BifurcationKit.AbstractContinuationAlgorithm,_contParams::ContinuationPar;kwargs...)
```

```@docs
BifurcationKit.continuation(prob_vf,
			bt::BifurcationKit.BogdanovTakens,
			bvp::BifurcationKit.AbstractBoundaryValueProblem,
			alg::BifurcationKit.AbstractContinuationAlgorithm,
			_contParams::ContinuationPar ;
			ϵ0 = 1e-5, amplitude = 1e-3,
			freeparams = ((@lens _.ϵ0), (@lens _.T)),
			maxT = Inf,
			update_every_step = 1,
			test_orbit_flip = false,
			test_inclination_flip = false,
			kwargs...
			)
```

## Utils 

```@docs
generate_hom_problem
```

```@docs
get_homoclinic_orbit
```

## Misc.
