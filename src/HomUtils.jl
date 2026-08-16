# change for Makie
function modify_hom_plot(::Union{BK.BK_NoPlot, BK.BK_Plots}, probHom, pars, lens; kwargs...)
    _plotsol = get(kwargs, :plot_solution, nothing)
    _plotsol2 = isnothing(_plotsol) ? BK.plot_default : (x, p; k...) -> _plotsol(x, (prob = probHom, lens = lens, p = p); k...)
end

function modify_hom_plot(::BK.BK_Makie, probHom, pars, lens; kwargs...)
    _plotsol = get(kwargs, :plot_solution, nothing)
    _plotsol2 = isnothing(_plotsol) ? BK.plot_default : (ax, x, p; k...) -> _plotsol(ax, x, (prob = probHom, lens = lens, p = p); k...)
end

# function to extract trajectories from branch
function get_homoclinic_orbit(br::BK.AbstractBranchResult, ind::Int)
    𝐇𝐨𝐦 = br.prob.VF.F
    x = br.sol[ind].x
    par0 = BK.setparam(br, br.sol[ind].p)
    get_homoclinic_orbit(𝐇𝐨𝐦, x, par0)
end