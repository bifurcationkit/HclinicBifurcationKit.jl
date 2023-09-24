# change for Makie
function modify_hom_plot(probPO, lens, kwargs)
    _plotsol = get(kwargs, :plot_solution, nothing)
    _plotsol2 = isnothing(_plotsol) ? (x, p; k...) -> nothing : (x, p; k...) -> _plotsol(x, (prob = probPO, lens = lens, p = p); k...)
end

# function to extract trajectories from branch
function get_homoclinic_orbit(br::BK.AbstractBranchResult, ind::Int)
    𝐇𝐨𝐦 = br.prob.VF.F
    x = br.sol[ind].x
    par0 = BK.setparam(br, br.sol[ind].p)
    get_homoclinic_orbit(𝐇𝐨𝐦, x, par0)
end