score <- function(models, D, control, verbose=TRUE, graphClass="graphNEL") {

  #if single model as input
  if (class(models)=="matrix") models <- list(models)    

  # Which Sgenes were silenced?
  Sgenes <- setdiff(unique(colnames(D)),"time")
  nrS <- length(Sgenes)
  
  # check that all models have S-genes as names
  fkt <- function(x,s){
     ss <- sort(s)
     c1 <- all(sort(colnames(x))==ss)
     c2 <- all(sort(rownames(x))==ss)
     return(c1 & c2)
  }  
  if (!all(sapply(models,fkt,s=Sgenes))) stop("\nnem:score> models must have same names as data")
  nrS <- length(Sgenes)     
  # make probability/density matrices D0 and D1  
  # nrow=#E-genes and ncol=#S-genes        
  if(control$type %in% c("mLL", "FULLmLL")){
  	D1 = sapply(Sgenes, function(s) rowSums(D[,colnames(D) == s,drop=FALSE]))  
  	D0 = sapply(Sgenes, function(s) sum(colnames(D) == s)) - D1 
  }
  else{
	D1 = D
	D0 = NULL
  }
    # if no prior is supplied:
  # assume uniform prior over E-gene positions      
  if (is.null(control$Pe)){ 	
	control$Pe <- matrix(1/nrS,nrow=nrow(D1),ncol=nrS)
	colnames(control$Pe) <- Sgenes  		
  }          
  if(control$type %in% c("CONTmLLRatio", "CONTmLLMAP")){		
  	control$Pe = cbind(control$Pe, double(nrow(D1)))  			
  	control$Pe[,ncol(control$Pe)] = control$delta/nrS
	control$Pe = control$Pe/rowSums(control$Pe)		
  }  
  if(is.null(control$Pm) & (control$lambda != 0)){
	cat(">>> Regularization parameter non-zero: Generating sparsity prior automatically! <<<\n")
	control$Pm = diag(length(Sgenes))
  }
   
    
  if (control$type=="FULLmLL"){ # FULL log marginal likelihood of all models
    if (verbose==TRUE) cat("Computing FULL (marginal) likelihood for",length(models),"models\n")
   
    results <- sapply(models,FULLmLL,D1,D0,control, verbose)        
  }
  else{   # log marginal likelihood of all models		
	if (verbose==TRUE) cat("Computing (marginal) likelihood for",length(models),"models\n")      		           
	results <- sapply(models,mLL,D1,D0,control, verbose)     	   			
  }     
  s       <- unlist(results["mLL",])
  ep      <- results["pos",]
  map     <- results["mappos",] 	
  LLperGene = results["LLperGene",]    
  para = results["para",]
#   if(!is.null(Pm)){  	
#   	log_pD_cond_Phi <- s  	  	  		
#   	if(is.null(control$lambda) || (control$lambda == 0)){  			
# 		if(verbose) cat("--> Using Bayesian model averaging to incorporate prior knowledge\n")
#   		lpPhi <- sapply(models, PhiDistr, Pm, a=1, b=0.5)		  		
#   		s <- log_pD_cond_Phi + lpPhi 		
#   	}
#   	else
# 		s = s - control$lambda*sapply(models, function(M) sum(abs(M - control$Pm))) + nrS^2*log(control$lambda*0.5)  
#   }    
  ppost = exp(s - max(s))
  ppost <- ppost/sum(ppost)  	# posterior model probability
  
  if(verbose){
	if(length(ppost) > 1){
		sort_ppost <- sort(ppost,decreasing=TRUE)		
		cat("(probabilities of best and second best model for ",Sgenes, ":", sort_ppost[1],",",sort_ppost[2],")\n")
	}
  }  
   # winning model       
  winner <- models[[which.max(s)]]  
  diag(winner) <- 0  
  if(graphClass == "graphNEL"){
  	gR <- new("graphAM",adjMat=winner,edgemode="directed")  
  	gR <- as(gR,"graphNEL")    
  }
  else
	gR <- winner
  PHI <- matrix(0,ncol=nrS,nrow=nrS)
  dimnames(PHI) <- list(Sgenes,Sgenes)  
  for(i in 1:length(models))
	PHI <- PHI + models[[i]]*ppost[i]   
  selected = results["mappos",which.max(s)][[1]]  
  selected = unique(unlist(selected[Sgenes]))  
  # output  
  res <- list(graph=gR, mLL=s, ppost=ppost, avg=PHI, pos=ep, mappos=map, control=control, selected=selected, LLperGene=LLperGene, para=para)
  class(res) <- "score"   
  return(res)  
}
