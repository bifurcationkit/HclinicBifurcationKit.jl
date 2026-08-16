# [Branch switching](@id Branch-switching-page)

```@contents
Pages = ["branchswitching.md"]
Depth = 3
```

## Branch switching from a Bogdanov–Takens (BT) point to a homoclinic curve

We provide an automatic branch switching method in this case (see for example the tutorial [Autonomous electronic circuit (aBS from BT)](@ref)). Hence, you can perform automatic branch switching by calling `continuation` with the following arguments:

```julia
continuation(prob_vf,
	bt,    # Bogdanov–Takens point
	bvp,   # discretization, e.g. Collocation or Shooting
	alg,   # continuation algorithm, e.g. PALC
	_contParams;
	ϵ0 = 1e-5, amplitude = 1e-3,
	maxT = Inf,
	kwargs...)
```

where `prob_vf` is the `BifurcationProblem` encoding the vector field and `bt` the Bogdanov–Takens point detected during a Fold or Hopf continuation, e.g. `bt = get_normal_form(br, ind_bt)`. More information is available in the docstring of [`continuation`](@ref) and in the [tutorials](@ref tutorials-page).
