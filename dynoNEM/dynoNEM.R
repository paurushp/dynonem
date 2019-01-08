# TODO: Add comment
# 
# Author: frohlich
###############################################################################
library(multicore) 
library(nem)
library(boot)
library(binom)

dyn.load("~/workingAt/trunk/dynoNEM/dynoNEMWrapper.so")

# simple DNEM-model, where time steps are assumed to be independent
# works only without replicates for discrete case!
simpleDNEM = function(D, priorNet=NULL, lambda=10^(-3:3), type=c("pval.dens", "effect.prob", "discrete"),  alpha=0.1, beta=0.2){
	T = dim(D)[1]
	nsgenes = dim(D)[3]
	negenes = dim(D)[2]
	if(is.null(priorNet))
		priorNet = diag(nsgenes)
	else{
		if(NCOL(priorNet) != nsgenes | NROW(priorNet) != nsgenes)
			stop("priorNet has to be a # S-genes x # S-genes matrix")
		if(any(priorNet > 1) | any(priorNet < 0))
			stop("All entries in priorNet have to be probabilities")
	}	
	if(any(lambda < 0))
		stop("Parameter lambda has to be >= 0")
	
	if(type == "pval.dens")
		D = log(D)
	if(type == "effect.prob")
		D = log(D/(1-D))
	myD = matrix(0, nrow=dim(D)[2], ncol=dim(D)[3])
	T = dim(D)[1]
	for(t in 1:T){		
		myD = myD + D[t,,]		
	}
	colnames(myD) = colnames(D[1,,])
	rownames(myD) = as.character(1:nrow(myD))
	if(!discrete)
		control = set.default.parameters(colnames(myD), type="CONTmLLMAP", Pm=priorNet, trans.close=FALSE)
	else
		control = set.default.parameters(colnames(myD), Pm=priorNet, para=c(alpha, beta))
	net = nemModelSelection(lambda, myD, control=control, verbose=FALSE)
	graph = as(net$graph, "matrix")
	diag(graph) = 0	
	list(network=graph, loglik=net$mLL)
}


dynoNEM = function(D, initial=NULL, priorNet=NULL, priorE=NULL, nu=10^(-3:3), type=c("pval.dens", "effect.prob", "discrete"), nrep=4, alpha=0.1, beta=0.2, bootstrap=FALSE, nboot=1000, mc.cores=1){	
	T = dim(D)[1]
	nsgenes = dim(D)[3]
	negenes = dim(D)[2]
	type = match.arg(type, several.ok=FALSE)
	if(type == "pval.dens")
		int.type = 0
	else if(type == "effect.prob")
		int.type = 1
	else if(type == "discrete")
		int.type = 2
	
	if(is.null(priorNet))
		priorNet = matrix(0, ncol=nsgenes, nrow=nsgenes)
	else{
		if(NCOL(priorNet) != nsgenes | NROW(priorNet) != nsgenes)
			stop("priorNet has to be a # S-genes x # S-genes matrix")
		if(any(priorNet > 1) | any(priorNet < 0))
			stop("All entries in priorNet have to be probabilities")
	}
	if(is.null(priorE))
		priorE = matrix(1/nsgenes, ncol=nsgenes, nrow=negenes)
	else{
		if(NCOL(priorE) != nsgenes | NROW(priorE) != negenes)
			stop("priorE has to be a # E-genes x # S-genes matrix")
		if(any(priorE > 1) | any(priorE < 0))
			stop("All entries in priorE have to be probabilities")
	}
	if(!is.null(initial)){
		if(NCOL(initial) != nsgenes | NROW(initial) != nsgenes)
			stop("initial has to be a # S-genes x # S-genes matrix")
		if(any(initial < 0))
			stop("All entries in initial have to be positive integers")
	}
	if(any(nu <= 0))
		stop("Parameter nu has to be positive")
	
	infer = function(idx.orig=1:negenes, idx){		
		Db = D[,idx,]		
		net = matrix(0, ncol=nsgenes, nrow=nsgenes)
		bics = mclapply(nu, function(n){								
					if(is.null(initial)){						
						if(type == "pval.dens")
							logD = log(Db)
						else if(type == "effect.prob")
							logD = log(Db/(1-Db))
						else
							logD = Db
						myD = matrix(0, nrow=dim(Db)[2], ncol=dim(Db)[3])					
						for(t in 1:T){
							myD = myD + logD[t,,]
						}						
						colnames(myD) = colnames(Db[1,,])
						rownames(myD) = as.character(1:nrow(myD))					
						if(type != "discrete")
							control = set.default.parameters(unique(colnames(myD)), type="CONTmLLMAP", Pm=priorNet, Pe=priorE, trans.close=FALSE, lambda=1/n)					
						else if(type == "discrete")
							control = set.default.parameters(unique(colnames(myD)), Pm=priorNet, Pe=priorE, lambda=1/n, para=c(alpha,beta))
						initial = nem(myD, control=control, verbose=FALSE)								
						initial = as(initial$graph, "matrix")
					}
					res = .C("dynoNEMWrapper", as.integer(dim(Db)[1]), as.integer(dim(Db)[3]), as.integer(dim(Db)[2]), D_R=Db, initial_R=initial, network_prior_R=priorNet, Egene_prior_R=priorE, n, as.integer(int.type), as.integer(nrep), alpha, beta, res_network=net, res_likelihood=double(1))
					BIC=-2*res$res_likelihood + log(nrow(Db[1,,]))*sum(abs(res$res_network - priorNet) > 0)				
					list(BIC=BIC, net=res$res_network, loglik=res$res_likelihood)					
				}, mc.cores=mc.cores)
		best = which.min(sapply(bics, function(bi) bi$BIC))
		net = bics[[best]]$net
		dimnames(net) = list(dimnames(D)[[3]], dimnames(D)[[3]])
		loglik = bics[[best]]$loglik
		if(!bootstrap)
			return(list(network=net, loglik=loglik, BIC=bics[[best]]$BIC, nu=nu[best]))
		else
			return(c(as.vector(net), loglik, bics[[best]]$BIC, nu[best]))
	}
	
	if(!bootstrap){
		net = infer(idx=1:negenes)
		return(net)
	}	
	else{
		set.seed(1234)
		net.boot = boot(1:negenes, infer, nboot)
		CIs = sapply((ncol(net.boot$t) - 2):ncol(net.boot$t), function(i) boot.ci(net.boot, type="normal", index=i))			
		sgenes = dimnames(D)[[3]]		
		med.times = apply(net.boot$t[,1:nsgenes^2], 2, median)
		med.times = matrix(med.times, ncol=nsgenes)
		dimnames(med.times) = list(sgenes, sgenes)
		CIs.net = binom.confint(colSums(net.boot$t[,1:nsgenes^2] > 0),nboot, methods="exact")[,c("lower", "upper")]
		CIs.net.lower = matrix(CIs.net[,"lower"], ncol=nsgenes)
		dimnames(CIs.net.lower) = list(sgenes, sgenes)
		CIs.net.upper = matrix(CIs.net[,"upper"], ncol=nsgenes)
		dimnames(CIs.net.upper) = list(sgenes, sgenes)			
		return(list(CIs.net.lower=CIs.net.lower, CIs.net.upper=CIs.net.upper, CIs=CIs, boot.result=net.boot, med.times=med.times))
	}	
}

dynoNEM.perturb = function(net, T, k){
	if(NROW(net) != NCOL(net))
		stop("net has to be a quadratic, weighted adjacency matrix")
	if(k > NROW(net) | k < 1)
		stop("k has to be between 1 and #S-genes")
	if(T < 0)
		stop("T has to be > 0")
	ret = .C("dynoNEM_getPerturbProb", net, as.integer(T), as.integer(NROW(net)), as.integer(k-1), res=matrix(0, ncol=T, nrow=nrow(net)))$res	
	rownames(ret) = rownames(net)
	colnames(ret) = as.character(1:T)
	ret
}

