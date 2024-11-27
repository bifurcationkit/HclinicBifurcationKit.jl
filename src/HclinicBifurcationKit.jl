module HclinicBifurcationKit
    using BifurcationKit, DocStringExtensions, Parameters
    using ForwardDiff
    using LinearAlgebra: norm, schur, ordschur, I, mul!, dot, eigen, normalize!
    using RecursiveArrayTools: ArrayPartition
    const BK = BifurcationKit



    include("HomProblemPBC.jl")
    include("HomUtils.jl")
    include("contKind.jl")

    include("PeriodicOrbitCollocation.jl")
    include("StandardShooting.jl")

    export generate_hom_problem, get_homoclinic_orbit

    BK.get_lenses(br::BK.AbstractResult{ <: HomoclinicHyperbolicSaddleCont}) = (getlens(br.prob), getlens(br.prob.VF.F))

end
