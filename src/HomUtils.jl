function modifyHomPlot(probPO, lens, kwargs)
	_plotsol = get(kwargs, :plotSolution, nothing)
	_plotsol2 = isnothing(_plotsol) ? (x, p; k...) -> nothing : (x, p; k...) -> _plotsol(x, (prob = probPO, lens = lens, p = p); k...)
end
