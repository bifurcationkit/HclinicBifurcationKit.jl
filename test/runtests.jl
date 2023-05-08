using HclinicBifurcationKit
using Test

@testset "HclinicBifurcationKit.jl" begin
	include("testHomoclinic.jl")
end

@testset "generate problem" begin
	include("generateHom.jl")
end
