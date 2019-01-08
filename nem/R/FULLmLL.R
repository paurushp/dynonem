# FULL log marginal likelihood of a models

FULLmLL <- function(Phi,D1,D0,control,verbose=FALSE) {
  if (!all(diag(Phi)==1)) diag(Phi) = 1
  
  # make sure model is in right order wrt to count matrices
  Phi <- Phi[colnames(D1),colnames(D1)]

  # compute score
  a0 <- control$hyperpara[1]
  b0 <- control$hyperpara[2]
  a1 <- control$hyperpara[3]
  b1 <- control$hyperpara[4]    
  n01 <- D1 %*% (1-Phi)
  n00 <- D0 %*% (1-Phi)
  n11 <- D1 %*% Phi
  n10 <- D0 %*% Phi
  s0  <- gamma(a0+b0)*gamma(n10+a0)*gamma(n00+b0)/gamma(a0)/gamma(b0)/gamma(n10+n00+a0+b0)
  s1  <- gamma(a1+b1)*gamma(n11+a1)*gamma(n01+b1)/gamma(a1)/gamma(b1)/gamma(n11+n01+a1+b1)  
  SP  <- s0*s1*control$Pe
  LLperGene = log(rowSums(SP))
  s   <- sum(LLperGene)
 
  # posterior effect positions
  ep  <- SP/rowSums(SP)

  # MAP estimate of effect positions
  Theta = apply(ep,1,function(e) e ==max(e))
  map = apply(Theta,1,which)	

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

  return(list(mLL=s,pos=ep,mappos=map,LLperGene=LLperGene, para=NULL))
}
