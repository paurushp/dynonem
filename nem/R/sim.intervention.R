sim.intervention = function(x, int){
	M = as(x$graph, "matrix")
	diag(M) = 1
	Sgenes.effected = lapply(int, function(i) colnames(M)[M[i,] == 1])
	Egenes.effected = x$mappos[unlist(Sgenes.effected)]
	list(Sgenes.effected=Sgenes.effected, Egenes.effected=Egenes.effected)
}
