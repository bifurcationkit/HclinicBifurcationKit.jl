module HclinicBifurcationKit
	using BifurcationKit, DocStringExtensions, Setfield, Parameters
	using ForwardDiff
	using LinearAlgebra: norm, schur, ordschur, I, mul!, dot, eigen, normalize!
	using RecursiveArrayTools: ArrayPartition
	const BK = BifurcationKit



	include("HomProblemPBC.jl")
	include("HomUtils.jl")
	include("contKind.jl")

	include("PeriodicOrbitCollocation.jl")
	include("StandardShooting.jl")

	export generateHomProblem, getHomoclinicOrbit

end