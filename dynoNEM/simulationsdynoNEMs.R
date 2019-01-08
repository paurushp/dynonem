library(RBGL)
library(Hmisc)
library(multicore)
library(nem)
source("~/workingAt/trunk/dynoNEM/dynoNEM.R")


## retrieve graphs of KEGG signaling pathways
KEGG.graphs = function(){	
	require(KEGGgraph)
	require(RCurl)
	require(gene2pathway)
	
	kegg_hierarchy = gene2pathway:::getKEGGHierarchy()
	signaltrans = kegg_hierarchy$pathIDsLev2[grep("Signal Transduction", names(kegg_hierarchy$pathIDsLev2))]
	pathways = names(which(sapply(kegg_hierarchy$parentPaths, function(p) signaltrans %in% unlist(p))))
	pathways = paste("hsa", pathways, ".xml", sep="")
	
	baseurl = "ftp://ftp.genome.jp/pub/kegg/xml/kgml/non-metabolic/organisms/hsa/"
	ftp = getURL(baseurl)
	ftp = unlist(strsplit(ftp,"\n"))
	ftp = ftp[grep("hsa[0-9][0-9][0-9][0-9][0-9]\\.xml", ftp)]
	urls = sapply(ftp, function(ft) paste(baseurl, substr(ft, nchar(ft)-11,nchar(ft)), sep=""))
	idx = unlist(sapply(pathways, function(p) grep(p, urls)))
	graphs = sapply(urls[idx], function(u) parseKGML2Graph(u))		
	nacc = sapply(graphs, function(mygraph) sapply(acc(mygraph, nodes(mygraph)), length))		
	save(graphs, nacc, file="~/workingAt/trunk/dynoNEM/KEGGgraphs.rda")
}

## DNEM model by Anchang et al. - very slow
#DNEM = function(D, initial=NULL, priorNet=NULL, priorE=NULL){
#	T = dim(D)[1]	
#	control = set.default.parameters(unique(colnames(D[T,,])), para=c(0.1, 0.1))	
#	if(!is.null(priorNet)){
#		control$Pm = priorNet
#		net = nemModelSelection(c(0.01, 0.01, 0.1, 1, 100, 100), D[T,,],inference="triples", control=control)
#	}
#	else
#		net = nem(D[T,,],inference="triples", control=control)
#	G = as(net$graph, "matrix")
#	D0<-array(0,dim=dim(D))
#	colnames(D0)<-colnames(D)
#	rownames(D0)<-rownames(D)
#	Data1<-abind(D0[1,,],D,along=1) 
#	theta.init = sapply(as.character(1:dim(D)[2]), function(i) sapply(net$mappos, function(g) i %in% g))
#	theta.init = apply(theta.init, 2, which.max)
#	lags.init=sample(0:T, dim(D)[2], replace=TRUE)
#	Data1 = aperm(Data1, c(2,3,1))
#	res = dnem(3, Data1, 10000, G, theta.init.lags.init, 0.1, 0.1, 0.2, file="./output.rda")
#	np = length(grep("lag", colnames(res[,,1])))
#	dnemoutput(G, file1="./output.rda", file2="./outputfinal.rda", burnin=1000, np=np,T=T, 0.4)	
#}

simulatePerturbation = function(Psi, T, k){	
	dynoNEM.perturb(Psi, T, k)
}

simulateData = function(net, T, decay=0.5, discrete=FALSE, nrepl=1){
    n = nrow(net$network)
    m = length(net$posEgenes)	
	D = array(0, dim=c(T, m, n))	
	for(k in 1:n){
		pertS = simulatePerturbation(net$network, T, k)
		for(t in 1:T){									
			effected = unique(which(net$posEgenes %in% which(pertS[,t]==1)))        
			not_effected = setdiff(1:m, effected)   
			if(!discrete){			
				b = sample(5:50,1)  # sample a random a
				a = sample(seq(0.1,0.9,by=0.1),1)
				l = c(sample(seq(0.01,0.49,by=0.01),1),sample(seq(0.01,0.49,by=0.01),1)) # sample mixing coefficient
				best = max(2,rnorm(1,b,b/10)) # put some noise on b
				aest = min(max(0.1,rnorm(1,a,0.05)),0.9) # and on a as well
				lest = pmax(0.01,l+rnorm(2,0,0.05))
				
				#cat("sample parameters: (alpha=",asamp,"beta=",bsamp,"lambda=",lsamp,")\n")					
				D[t, effected,k] =  nem:::bum.ralt(length(effected),c(a,b),c(1-sum(l),l[1],l[2])) #sample p-values for effected genes => aus H1 ziehen
				D[t, not_effected,k] = runif(length(not_effected))  # ... and not effected ones => aus H0 ziehen   			
				D[t,,k] = nem:::bum.dalt(D[t,,k],c(aest,best),c(1-sum(lest),lest[1],lest[2]))								
			}
			else{
				a = sample(c(0.01, 0.05, 0.1, 0.2, 0.3), 1)
				b = sample(c(0.01, 0.05, 0.1, 0.2, 0.3), 1)		
				D[t, effected, k] = rbinom(length(effected), nrepl, p=(1-a))
				D[t, not_effected, k] = rbinom(length(not_effected), nrepl, p=b)
			}				
		}
	}
	if(!discrete)
		D[D == 0] = min(D[D > 0]) 
	sgenes = paste("S",seq(1,n,1),sep="")
	dimnames(D) = list(as.character(1:T), as.character(1:m), sgenes)
    D
}

sampleKEGGPathway = function(graphs, nacc, n){
	mygraph = sample(graphs, 1)
	mynacc = nacc[[names(mygraph)]]
	mygraph = mygraph[[1]]
	u_mygraph = ugraph(mygraph)		
	start = names(sample(mynacc[mynacc >= n],1))
	mynodes = start
	no = start
	while(length(mynodes) < n){
		nei = adj(u_mygraph, no)
		no = sample(nei[[1]], 1)[[1]]           
		if(!(no %in% mynodes))
			mynodes = c(mynodes, no)
	}
	S = subGraph(mynodes, mygraph)		
	S = removeSelfLoops(S)	
	S
}

sampleRandomGraph = function(n, p=2/n){
	require(igraph)
	S = erdos.renyi.game(n, p, directed=T)
	S = get.adjacency(S)
	S = as(S, "graphNEL")
	S
}

sample.DNEM.structure = function(S, m, Tmax=4, attach="uniform", decay=0.5, prob=0.5){	
	n = ncol(S)	
	sgenes = paste("S",seq(1,n,1),sep="")
	colnames(S) = sgenes
	rownames(S) = sgenes
	Strans = nem:::transitive.closure(S, mat=T, loops=F)		
	Strans[Strans==1] = pmin(rgeom(sum(Strans), prob=prob) + 1, Tmax)	
	
	# correct time lags: indirect edges have to be faster than direct ones in order to be observable
	for(y in 1:nrow(Strans)){
		for(x in 1:nrow(Strans)){
			if(Strans[x,y] != 0){
				for(j in 1:nrow(S)){
					if((Strans[y,j] != 0) && (Strans[x,j] > 1) && (Strans[x,j] >= Strans[y,j] + Strans[x,y]))
						Strans[x,j] = Strans[y,j] + Strans[x,y] - 1						
				}
			}
		}
	}
	edges = which(S == 1)
	S[edges] = Strans[edges]		
	
    if(attach=="uniform") # uniform distribution of E-genes
        epos = sample(1:n,m,replace=TRUE) # position of E-genes uniform             
    else if(attach %in% c("downstream","upstream")){ # E-genes preferentially downstream or upstream        
        g = as(S,"graphNEL")        
        startnode = which.min(apply(S,2,sum)) # start with node with minimal in-degree
        visit = bfs(g,nodes(g)[startnode])      
        dis = 0:(n-1)       
        names(dis) = visit
        r = 1 - decay*dis/(n-1)
        r = r/sum(r)
        if(attach == "downstream")
            r = rev(r)
        epos = sample(1:n,m,replace=TRUE,prob=r)
    }           
    list(network=S,posEgenes=epos)      
}

networkInference = function(D, net, method="dynoNEM", usePrior="sparse", fracErr=0, fracKnown=0.25, discrete=FALSE){ # perform network inference and evaluate the best model         
    print(usePrior)     
    original = (net$network > 0)*1
    posEgenes = net$posEgenes
	nrS = NCOL(original)
       
    if(usePrior != "no"){   
        if(length(grep("sparse", usePrior) > 0))
            priorPHI = diag(nrS)
        else{
            priorPHI = matrix(0.4,ncol=nrS,nrow=nrS)
            dimnames(priorPHI) = dimnames(original)     
            known = which(original == 1)
            nEdges = length(known)
            r = sample(1:nEdges, floor(fracKnown*nEdges))       
            priorPHI[known[r]] = 1
            unknown = which(original == 0)                        
            r = sample(1:length(unknown), floor(fracErr*length(unknown)))               
            priorPHI[unknown[r]] = 1                            
            print((priorPHI>0.5)*1)
            cat("total known edges = ",sum(priorPHI==1),"\n")           
        }
             
    }
    else
        priorPHI = NULL             
       				
        			
	if(method == "dynoNEM")
   		elapsed=system.time(for(i in 1:1) est = dynoNEM(D, nu=c(0.01, 0.1, 1, 10, 100), discrete=discrete, nrep=1))[1]
	else if(method =="simpleDNEM")
		elapsed=system.time(for(i in 1:1) est = simpleDNEM(D, discrete=discrete))[1] # simple nem model
	else if(method == "DNEM")
		elapsed=system.time(for(i in 1:1) est = DNEM(D, priorNet=priorPHI))[1] # DNEM	
	else
		stop("unknown method")
   	inferred = abs(est$network)
	mse = mean((inferred - net$network)^2)
   	mat = (inferred > 0)*1
               
    print("original:")      
    print(net$network)
    print("inferred:")    
    print(inferred)          
    tp = sum(mat[original==1]==1,na.rm=TRUE)
    tn = sum(mat[original==0]==0,na.rm=TRUE)
	if(discrete & method == "simpleDNEM")
		fp = sum(mat[transitive.closure(original, mat=TRUE, loops=F)==0]==1,na.rm=TRUE)
	else
		fp = sum(mat[original==0]==1,na.rm=TRUE)
    fn = sum(mat[original==1]==0,na.rm=TRUE)    
    tpr = tp/(tp+fn)
    fpr = fp/(fp+tn)	
    if(is.nan(tpr))
        tpr = 1 		
    cat("tpr = ",tpr,"\n")
    cat("fpr = ",fpr,"\n")    
	cat("mse = ", mse, "\n")
    print(paste("elapsed time (s) = ", elapsed))    
    list(tpr=tpr,fpr=fpr,mse=mse)
}

# x: either # E-genes or number of time points for a fixed number of E-genes (4*n)
# usePrior: which prior to use ("no" means pure ML estimate)
# xlabel: put '# time points', if dependency on # time points is computed
# fracKnown: number of known edges in prior network (not needed here)
test = function(n, x=c(1, 2, 5, 10, 20)*n, usePrior="sparse", xlabel="# E-genes", fracKnown=0, methods=c("dynoNEM","simpleDNEM"), discrete=FALSE, outputdir="."){     
	set.seed(123456789)
	load("~/workingAt/trunk/dynoNEM/KEGGgraphs.rda")
	pdf(file.path(outputdir, paste("sampledNetworks_",n,"Sgenes.pdf")))
	subnets = list()
	i = 0
	while(i < 10){
		cand.net = sampleKEGGPathway(graphs, nacc, n)
		cand.net.mat = as(cand.net, "matrix")
		exists.net = any(sapply(subnets, function(subn) all(subn == cand.net.mat))) # CAUTION: This is not an exact test for graph isomorphism!!!
		if(!exists.net && (i < 9 || (i == 9 && length(tsort(cand.net)) == 0))){ # at least one graph has to have a cycle			
			i = i + 1			
			plot(cand.net)
			subnets[[i]] = cand.net.mat			
		}
	}
	dev.off()	
	ntrials = 100
    if(usePrior == "no" | usePrior=="sparse")
        fracKnown = 0         
	results = array(0, dim=c(10, 3, ntrials, length(x), length(methods)))	
	dimnames(results)[[2]] = c("tpr", "fpr", "mse")
	dimnames(results)[[5]] = methods		
	for(s in 1:length(subnets)){
		g = subnets[[s]]
		for(m in 1:length(x)){              
			res = mclapply(1:ntrials, function(i){
				restmp = array(dim=c(3, length(methods)))    	
				if(xlabel == "# E-genes"){
	            	net = sample.DNEM.structure(g, x[m])
					D = simulateData(net, T=10, discrete=discrete)
				}
				else if(xlabel == "# time points"){
					net = sample.DNEM.structure(g, 10*n)			
	            	D = simulateData(net, T=x[m], discrete=discrete)
				}
				else if(xlabel == "parameter p"){
					net = sample.DNEM.structure(g, 10*n, prob=x[m])			
					D = simulateData(net, T=10, discrete=discrete)
				}
				else
					stop("unknown test procedure")
				for(method in 1:length(methods)){
	            	res = networkInference(D, net, methods[method],usePrior=usePrior, fracKnown=fracKnown, discrete=discrete)              
	            	restmp[,method] = unlist(res)  
				}
				restmp
	        }, mc.cores=7) 
			for(i in 1:ntrials)
				results[s,,i,m,] = res[[i]] 
	
			for(method in 1:length(methods)){
				cat("method: ", methods[method],"\n\n")
				print(paste("--> mean % tpr (", xlabel, " =", x[m],") = ", rowMeans(results[s,,,m,method])[1]))           
				print(paste("--> mean % fpr (", xlabel, " =", x[m],") = ", rowMeans(results[s,,,m,method])[2]))
				print(paste("--> mean MSE (", xlabel, " =", x[m],") = ", rowMeans(results[s,,,m,method])[3]))
				cat("====================================\n")
			}		
	    }
	}
    save(results,file=file.path(outputdir, paste("results_n",n, "_", xlabel, "_", usePrior, "_", fracKnown, ".rda",sep="")))
	
    pdf(file.path(outputdir, paste("sensitivity_n",n, "_", xlabel, "_", usePrior, "_", fracKnown, ".pdf", sep="")))    
	plot(x,seq(0,100,length.out=length(x)),type="n",main=paste("n = ",n),xlab=xlabel, ylab="sensitivity (%)")
	for(method in 1:length(methods)){
		means = 100*apply(apply(results[,1,,,method],c(3,1), mean), 1, mean)
		sds = 100*apply(apply(results[,1,,,method],c(3,1), mean), 1, sd)
		lines(x, means, lty=method, lwd=2, type="b")
    	errbar(x,means,means+sds,means-sds,xlab=xlabel,ylab="sensitivity (%)",main=paste("n = ",n),ylim=c(0,100),type="p", add=TRUE, lwd=2)
	}
	legend("bottomright", methods, lty=1:length(methods))
    dev.off()
	
	pdf(file.path(outputdir, paste("specificity_n",n, "_", xlabel, "_", usePrior, "_", fracKnown, ".pdf", sep="")))    
	plot(x,seq(0,100,length.out=length(x)),type="n",main=paste("n = ",n),xlab=xlabel, ylab="1 - fpr (%)")
	for(method in 1:length(methods)){		
		means = 100*(1 - apply(apply(results[,2,,,method],c(3,1), mean), 1, mean))
		sds = 100*apply(apply(results[,2,,,method],c(3,1), mean), 1, sd)		
		lines(x, means, lty=method, lwd=2, type="b")
		errbar(x,means,means+sds,means-sds,xlab=xlabel,ylab="1 - fpr (%)",main=paste("n = ",n),ylim=c(0,100),type="p", add=TRUE, lwd=2)
	}
	legend("bottomright", methods, lty=1:length(methods))
	dev.off()
	
	pdf(file.path(outputdir, paste("MSE_n",n, "_", xlabel, "_", usePrior, "_", fracKnown, ".pdf", sep="")))    
	method = which(methods=="dynoNEM")
	plot(x,seq(0,max(results[,3,,,method]),length.out=length(x)),type="n",main=paste("n = ",n),xlab=xlabel, ylab="MSE")	
	#for(method in 1:length(methods)){
		means = apply(apply(results[,3,,,method],c(3,1), mean), 1, mean)
		sds = apply(apply(results[,3,,,method],c(3,1), mean), 1, sd)
		lines(x, means, lty=method, lwd=2, type="b")
		errbar(x,means,means+sds,means-sds,xlab=xlabel,type="p", add=TRUE, lwd=2)
	#}
	#legend("bottomright", methods, lty=1:length(methods))
	dev.off()
		
	results
}

analyze.topologies = function(){
	require(ggplot2)	
	load("/home/bit/frohlich/workingAt/trunk/dynoNEM/simresults/results_n5_# time points_sparse_0.rda")
	times = c(3, 4, 6, 8, 10)
	for(T in 1:5){
		sens = data.frame(sensitivity=rbind(results[,1,,T,1], results[,1,,T,2]), method=as.factor(rep(c("dynoNEM", "simpleDNEM"), each=10)), network=as.factor(rep(1:10,2)))	
		sens = reshape(sens, varying=grep("sensitivity", colnames(sens)), times=1:100, timevar="trial", v.names="sensitivity", direction="long")
		spec = data.frame(specificity=1-rbind(results[,2,,T,1], results[,2,,T,2]), method=as.factor(rep(c("dynoNEM", "simpleDNEM"), each=10)), network=as.factor(rep(1:10,2)))	
		spec = reshape(spec, varying=grep("specificity", colnames(spec)), times=1:100, timevar="trial", v.names="specificity", direction="long")
		qplot(network, sensitivity, data=sens, geom="boxplot", fill=method) + theme_bw() + scale_fill_grey(end=1, start=0.5) + opts(axis.text.x = theme_text(hjust=1, angle=45), strip.text.x = theme_text(angle=90))
		ggsave(paste("/home/bit/frohlich/workingAt/trunk/dynoNEM/simresults/results_n5_networkArchitecture_sens_T",times[T],".pdf",sep=""))
		qplot(network, specificity, data=spec, geom="boxplot", fill=method) + theme_bw() + scale_fill_grey(end=1, start=0.5) + opts(axis.text.x = theme_text(hjust=1, angle=45), strip.text.x = theme_text(angle=90))
		ggsave(paste("/home/bit/frohlich/workingAt/trunk/dynoNEM/simresults/results_n5_networkArchitecture_spec_T",times[T],".pdf", sep=""))
	}
}
