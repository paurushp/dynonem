nemModelSelection <- function(lambdas,D,inference="nem.greedy",models=NULL,control=set.default.parameters(unique(colnames(D))),verbose=TRUE,...){
    infer <- function(lam){                     
        control$lambda=lam          
        res <- nem(D,inference=inference,models=models,control=control, verbose=verbose) # ACHTUNG: nem spuckt immer den MAP score aus!!!       
        if(length(res$selected) <= 1)
            res$mLL = -Inf
        else if(!is.null(control$Pm)){
            controltmp = control
            controltmp$Pm = NULL
            controltmp$lambda=0 
            if(control$type != "depn")
                resmLL <- score(list(as(res$graph,"matrix")),D[res$selected,,drop=FALSE],controltmp,verbose=FALSE)$mLL # get true mLL
            else
                resmLL <- score(list(as(res$graph,"matrix")),D,controltmp,verbose=FALSE)$mLL # get true mLL
            res$mLL = resmLL                
        }
        return(res)
    }   
    results <- lapply(lambdas,infer)            
    AICs <- lapply(results,network.AIC,verbose=verbose,control$Pm,k=log(nrow(D)),...)               
    winner <- results[[which.min(AICs)]]    
    if(verbose)
        cat(paste("====> chosen best model with lambda =",winner$control$lambda,"\n"))      
    return(winner)
}
