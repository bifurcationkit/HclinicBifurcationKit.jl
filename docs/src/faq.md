# FAQ

For general questions on continuation, bifurcation detection or the use of `BifurcationKit`, see also the [FAQ of BifurcationKit.jl](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/faq/).

## Common questions

- **How do I compute a homoclinic orbit from scratch?** You first need a periodic orbit (computed with BifurcationKit), then you build the homoclinic problem with `generate_hom_problem` and continue it; see the [tutorials](@ref tutorials-page).
- **How do I detect codimension-two bifurcations?** Enable the detection options during the continuation of the homoclinic branch; see [Detection of bifurcation points](@ref).
- **Which discretizations are supported?** Orthogonal collocation and (multiple) standard shooting, see [Homoclinic based on orthogonal collocation](@ref) and [Homoclinic based on parallel multiple shooting](@ref).

Do not hesitate to open an issue on [GitHub](https://github.com/bifurcationkit/HclinicBifurcationKit.jl/issues) if you have a question or a bug to report.
