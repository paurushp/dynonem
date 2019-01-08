# log marginal likelihood of model
mLL <- function(Phi,D1,D0=NULL,control, verbose=FALSE) {       
   if (!all(diag(Phi)==1)) 
  	diag(Phi) <- 1 	        
  para=NULL
  if(control$type == "mLL")
  	L  <- control$para[1]^(D1 %*% (1-Phi)) * (1-control$para[1])^(D0 %*% (1-Phi)) * (1-control$para[2])^(D1 %*% Phi) * control$para[2]^(D0 %*% Phi)  	    
  else if(control$type %in% c("CONTmLLBayes", "CONTmLLDens"))	
	L <- exp(D1%*%Phi)
  else if(control$type == "CONTmLL")
	L <- exp(log(D1)%*%Phi + log((1-D1))%*%(1-Phi)) 
  else if(control$type %in% c("CONTmLLMAP", "CONTmLLRatio")){			
	Phi2 = cbind(Phi, double(ncol(Phi)))
	colnames(Phi2)[ncol(Phi2)] = "null"
	ep = D1%*%Phi2 + log(control$Pe)	
	Theta = apply(ep,1,function(e) e ==max(e))
	L = t((Phi2%*%(Theta*1)>0)*1)*D1
	LLperGene=rowSums(L)		
	s = sum(LLperGene)			
	map = apply(Theta,1,which)	
  }
  else if(control$type == "depn"){	
	res = score.network(D1, Phi, control)
	s = res$loglik
	LLperGene = res$LLperSample
	ep = effect.likelihood(D1, res$net)		
	Theta = apply(ep,1,function(e) e ==max(e))
	map = apply(Theta,1,which)
	para = res$net$parameters
  }
  if(!(control$type %in% c("CONTmLLMAP","CONTmLLRatio", "depn"))){	
	if(!is.null(control$Pe))
		LP <- L*control$Pe
	else
		LP <- L	
	LLperGene = log(rowSums(LP))
	ep <- LP/(rowSums(LP))  	  			
	Theta = apply(ep,1,function(e) e ==max(e))
	s  <- sum(LLperGene)
  	map = apply(Theta,1,which)			
  }      
  if(!is.null(rownames(D1)))
	map = sapply(map, names)    
  if(!is.null(control$Pm)){
  	if(control$lambda != 0){
		if(verbose) cat("--> Using regularization to incorporate prior knowledge\n")			
  		s <- s - control$lambda*sum(abs(Phi - control$Pm)) + ncol(Phi)^2*log(control$lambda*0.5)    
	}
	else{
		if(verbose) cat("--> Using Bayesian model averaging to incorporate prior knowledge\n")
		s = s + PhiDistr(Phi, control$Pm, a=1, b=0.5)
	}
  }
  list(mLL=s,pos=ep,mappos=map,LLperGene=LLperGene, para=para)
}

