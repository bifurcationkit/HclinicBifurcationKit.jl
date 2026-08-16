using Pkg

using Documenter, HclinicBifurcationKit, BifurcationKit
ENV["GKSwstype"] = "100"

# to display progress
ENV["JULIA_DEBUG"] = Documenter

makedocs(doctest = false,
	sitename = "Homoclinic / Heteroclinic orbits in Julia",
	authors = "Romain Veltz",
	format = Documenter.HTML(
		collapselevel = 1,
		assets=[
			asset("https://bifurcationkit.github.io/assets/js/documentation.js"),
			asset("https://bifurcationkit.github.io/assets/css/documentation.css"),
				],
	),
	pages = Any[
		"Home" => "index.md",
		"Tutorials" => "tutorials/tutorials.md",
		"Problems" => [
			"Homoclinic to hyperbolic saddle" => [
				"Introduction" => "homoHS.md",
				"Collocation" => "periodicOrbitCollocation.md",
				"Shooting" => "shooting.md",
				],
		],
		"Functionalities" => [
			"Bifurcations" => [
				"Bifurcation detection (codim 1)" => "detectionBifurcation.md",
				"Branch switching" => "branchswitching.md",
							  ],
		],
		"Frequently Asked Questions" => "faq.md",
		"Library" => "library.md"
	]
	)

deploydocs(
	repo = "github.com/bifurcationkit/HclinicBifurcationKit.jl.git",
	devbranch = "main"
)
