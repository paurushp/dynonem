transitive.reduction <- function(g){
	if (!(class(g)%in%c("matrix","graphNEL"))) stop("Input must be an adjacency matrix or graphNEL object")
	if(class(g) == "graphNEL"){
		g = as(g, "matrix")		
	}
# 	if(class(g) == "matrix"){		
		# modified algorithm from Sedgewick book: just remove transitive edges instead of inserting them		
		g = g - diag(diag(g))		
		type = (g > 1)*1 - (g < 0)*1		
		for(y in 1:nrow(g)){
			for(x in 1:nrow(g)){
				if(g[x,y] != 0){
					for(j in 1:nrow(g)){
						if((g[y,j] != 0) && (g[x,j] != 0) & (sign(type[x,j])*sign(type[x,y])*sign(type[y,j]) != -1))
							g[x,j] = 0						
					}
				}
			}
		}
# 	}
# 	else{		
# 		nodenames=nodes(g)		
# 		for(y in 1:length(nodes(g))){
# 			edges = edgeL(g)
# 			x = which(sapply(edges,function(l) y %in% unlist(l)))
# 			j = unlist(edges[[y]])
# 			cands = sapply(edges[x], function(e) list(intersect(unlist(e),j)))			
# 			cands = cands[sapply(cands,length) > 0]
# 			if(length(cands) > 0)
# 				for(c in 1:length(cands)){ 				
# 					jj = unlist(cands[c])					
# 					g = removeEdge(rep(names(cands)[c],length(jj)),nodenames[jj],g)
# 				}
# 		}
# 	}	
	g		
}