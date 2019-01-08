set.default.parameters = function(Sgenes, ...){         
    map = as.list(Sgenes)
    names(map) = Sgenes 
    control = list(Sgenes=Sgenes,type="mLL",para=c(0.13,0.05),hyperpara=c(1,9,9,1),Pe=NULL,Pm=NULL,Pmlocal=NULL,local.prior.size=length(Sgenes),local.prior.bias=1,triples.thrsh=0.5,lambda=0,delta=1,selEGenes=FALSE, trans.close=TRUE, backward.elimination=FALSE, mode="continous_Bayesian", lambda.intervention=4, lambda.no_intervention=4,  df.intervention=4.4, df.no_intervention=4.4, nu.intervention=0.6, nu.no_intervention=0.95, scale.intervention=0.023, scale.no_intervention=0.023, map=map, outputdir="QualityControl", debug=FALSE)
    args = list(...)
    na = names(args)
    for(a in na){
        control[[a]] = args[[match(a, na)]]
    }
    control
}

nem <- function(D,inference="nem.greedy",models=NULL,control=set.default.parameters(setdiff(unique(colnames(D)),"time")), verbose=TRUE){
#------------------------------
# Sanity checks              
if(inference == "nem.greedyMAP")
    type = "CONTmLLMAP"
  
if (!(inference %in% c("pairwise", "triples", "search","ModuleNetwork","nem.greedy","nem.greedyMAP","BN.greedy","BN.exhaustive"))) 
    stop("\nnem> argument 'inference' is not valid\n")
if (!(control$type %in% c("mLL", "FULLmLL", "CONTmLL","CONTmLLBayes","CONTmLLMAP", "CONTmLLDens", "CONTmLLRatio","depn"))) 
    stop("\nnem> argument 'type' is not valid")
if (is.null(control$para) & is.null(control$hyperpara) & !(control$type %in% c("CONTmLL", "CONTmLLBayes", "CONTmLLMAP", "CONTmLLDens", "CONTmLLRatio","depn")))
    stop("\nnem> provide either 'para' or 'hyperpara'\n")
if (control$type == "mLL" & is.null(control$para)) 
    stop("\nnem> provide argument 'para'\n")
if (control$type == "FULLmLL" & is.null(control$hyperpara)) 
    stop("\nnem> provide argument 'hyperpara'\n")
if (!is.null(control$hyperpara)) {
    if (length(control$hyperpara) != 4) 
        stop("\nnem> 'hyperpara' is not a vector of length 4")
    if (!all(control$hyperpara > 0)) 
        stop("\nnem> 'hyperpara' must be >0")
}
if (!is.null(control$para)) {
    if (length(control$para) != 2) 
        stop("\nnem> 'para' is not a vector of length 2")
    if (any(control$para < 0) | any(control$para > 1)) 
        stop("\nnem> 'para' must be in [0,1]")
}
if(control$lambda < 0) control$lambda <- abs(control$lambda)

Sgenes <- unique(colnames(D))
if(control$selEGenes){  
    if(inference %in% c("BN.greedy","BN.exhaustive","depn"))
        stop("No automatic feature selection implemented for this inference method so far!\n")
    return(nem.featureselection(D, inference, models, control, verbose))
}

#------------------------------
# PAIRWISE                     

if (inference == "pairwise"){   
    if (is.null(control$local.prior.size))
        Pmlocal <- NULL 
    else{ 
        if (control$local.prior.size <= 0 | control$local.prior.bias <= 0) 
            stop("\nnem> local prior parameters invalid")
            Pmlocal <- local.model.prior(control$local.prior.size,length(Sgenes),control$local.prior.bias)      
        }
        result <- pairwise.posterior(D,control,verbose)    
    }

#------------------------------
# MODULE NETWORK                       
else if(inference == "ModuleNetwork"){
    result <- moduleNetwork(D,control,verbose=verbose)  
}

#------------------------------
# TRIPLES                     

else if (inference == "triples"){

        result <- triples.posterior(D,control,verbose)
        #A.t <- transitive.lp(A)

    }

#------------------------------
# GREEDY                     
else if(inference == "nem.greedy"){
    result <- nem.greedy(D,initial=models[[1]],control, verbose=verbose)    
}

else if(inference == "nem.greedyMAP"){
    result <- nem.greedyMAP(D,control, verbose=verbose) 
}

#------------------------------
# SEARCH                       

else if (inference == "search"){    
        if (is.null(models)) models <- enumerate.models(length(Sgenes),Sgenes,trans.close=control$trans.close,verbose=verbose)
        result <- score(models,D,control,verbose)
}

# Bayesian network
else if(inference %in% c("BN.greedy","BN.exhaustive")){
    if(inference == "BN.greedy")
        result = nem.BN(D, "greedy", mode=control$mode, Pm=control$Pm, lambda=control$lambda, verbose=verbose)
    else
        result = nem.BN(D, "exhaustive", mode=control$mode, Pm=control$Pm, lambda=control$lambda, verbose=verbose)
}
else
    stop(paste("Unknown inference method", inference,"\n"))

#------------------------------
# OUTPUT                       
#class(result) <- "nem"
return(result)

}
