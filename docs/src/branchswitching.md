# [Branch switching](@id Branch-switching-page)


```@contents
Pages = ["branchswitching.md"]
Depth = 3
```

## Branch switching from Bogdanov-Takens (BT) point to Homoclinic curve curve

We provide an automatic branch switching method in this case (see for example [Autonomous electronic circuit](@ref). Hence, you can perform automatic branch switching by calling `continuation` with the following options:

```julia
continuation(prob_vf, 
	# bt point
	bt,
	# also like ShootingProblem or PeriodicOrbitOCollProblem
	bvp, 
	# PALC, etc
	alg, 
	_contParams; 
	kwargs...)
```

Note that the BT has been detected during Fold or Hopf continuation. More information is available at [`continuation`](@ref)