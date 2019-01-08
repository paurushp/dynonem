nem.bootstrap <- function(D, thresh=0.5, nboot=1000,inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))), verbose=TRUE){
		
	inferNetwork <- function(D, boot){				
		controltmp = control				
		controltmp$Pe = control$Pe[boot,]
		Dtmp = D[boot,]
		if(!is.null(control$Pm) & length(control$lambda) > 1)
			res = as.vector(as(nemModelSelection(control$lambda,Dtmp,inference,models,controltmp,verbose)$graph,"matrix"))
		else
			res = as.vector(as(nem(Dtmp,inference,models,controltmp,verbose)$graph,"matrix"))
		res
	}
	#results = bootstrap(1:nrow(D),nboot,theta=inferNetwork)$thetastar
	if(!is.null(rownames(D))){
		if("time"  %in% colnames(D))
			group = as.factor(paste(rownames(D), D[,"time"],sep=""))
		else
			group = as.factor(rownames(D))
		res.boot = boot(D, inferNetwork, nboot, strata=group)
	}
	else
		res.boot = boot(D, inferNetwork, nboot)
	results =  res.boot$t
	Sgenes = setdiff(unique(colnames(D)), "time")
	n = length(Sgenes)
	overlapBoot = colMeans(results)
	overlapBoot = matrix(round(overlapBoot,digits=2),ncol=n,nrow=n)
	colnames(overlapBoot) = Sgenes
	rownames(overlapBoot) = Sgenes
	res = nem(D,models=list((overlapBoot>thresh)*1),inference="search",control, verbose=verbose)
	res$pos = res$pos[[1]]
	res$mappos = res$mappos[[1]]
	res$mLL = res$mLL[[1]]
	res$LLperGene = res$LLperGene[[1]]
	res$para = res$para[[1]]
	res$control= control
	g = res$graph
	edgeDataDefaults(g, "label") = 1	
	edgeDataDefaults(g, "weight") = 1
	for(s1 in Sgenes){
		for(s2 in Sgenes){
			if(s2 %in% unlist(adj(g, s1))){
				edgeData(g, from = s1, to = s2, attr = "weight") = overlapBoot[s1,s2]			
				edgeData(g, from = s1, to = s2, attr = "label") = overlapBoot[s1,s2]
			}
		}
	}
	res$graph = g
	class(res) <- "nem.bootstrap"
	res
}
