\name{nem}
\alias{nem}
\alias{print.nem}
\alias{print.nem.greedy}
\alias{print.nem.greedyMAP}
\alias{print.pairwise}
\alias{print.triples}
\alias{print.ModuleNetwork}
\alias{print.score}
\alias{print.nem.BN}

\title{Nested Effects Models - main function}
\description{The main function to perform model learning from data}
\usage{
nem(D,inference="nem.greedy",models=NULL,control=set.default.parameters(setdiff(unique(colnames(D)),"time")), verbose=TRUE)

\method{print}{nem}(x, ...)
}

\arguments{
  \item{D}{data matrix with experiments in the columns (binary or continious)}
  \item{inference}{\code{search} to use exhaustive enumeration, \code{triples} for triple-based inference, \code{pairwise} for the pairwise heuristic, \code{ModuleNetwork} for the module based inference, \code{nem.greedy} for greedy hillclimbing, \code{nem.greedyMAP} for alternating MAP optimization using log odds or log p-value densities}
  \item{models}{a list of adjacency matrices for model search. If NULL, an  exhaustive enumeration of all possible models is performed.}
  \item{control}{list of parameters: see \code{set.default.parameters}}
  \item{verbose}{do you want to see progression statements? Default: TRUE}

  \item{x}{nem object}
  \item{...}{other arguments to pass}
}
\details{  
   If parameter Pm != NULL and parameter lambda == 0, a Bayesian approach to include prior knowledge is used. Alternatively, the regularization parameter lambda can be tuned in a model selection step via the function \code{nemModelSelection} using the BIC criterion. 
   If automated subset selection of effect reporters is used and parameter type == CONTmLLMAP, the regularization parameter delta is tuned via the AIC model selection criterion. Otherwise, an iterative algorithm is executed, which in an alternating optimization scheme reconstructs a network given the current set of effect reporters and then selects the effect reporters having the highest likelihood under the given network. The procedure is run until convergence.

  The function \code{plot.nem} plots the inferred phenotypic hierarchy as a directed graph, the likelihood distribution of the models (only for exhaustive search) or the posterior position of the effected genes.  
}
\value{
  \item{graph}{the inferred directed graph (graphNEL object)}   
  \item{mLL}{log posterior marginal likelihood of final model} 
  \item{pos}{posterior over effect positions}
  \item{mappos}{MAP estimate of effect positions}   
  \item{selected}{selected E-gene subset}
  \item{LLperGene}{likelihood per selected E-gene}
  \item{control}{hyperparameter as in function call}
}

\author{Holger Froehlich <URL: http:/www.dkfz.de/mga2/people/froehlich>, Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}


\seealso{\code{\link{set.default.parameters}}, \code{\link{nemModelSelection}}, \code{\link{nem.jackknife}}, \code{\link{nem.bootstrap}}, \code{\link{nem.consensus}}, \code{\link{local.model.prior}}, \code{\link{plot.nem}}}
\examples{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]
   control = set.default.parameters(unique(colnames(D)), para=c(0.13, 0.05))   
   res1 <- nem(D,inference="search", control=control)
   res2 <- nem(D,inference="pairwise", control=control)
   res3 <- nem(D,inference="triples", control=control)
   res4 <- nem(D,inference="ModuleNetwork", control=control)
   res5 <- nem(D,inference="nem.greedy", control=control)        
   res6 = nem(BoutrosRNAiLods, inference="nem.greedyMAP", control=control)
   

   par(mfrow=c(2,3))
   plot.nem(res1,main="exhaustive search")
   plot.nem(res2,main="pairs")
   plot.nem(res3,main="triples")
   plot.nem(res4,main="module network")
   plot.nem(res5,main="greedy hillclimber")      
   plot.nem(res6,main="alternating MAP optimization")
}
\keyword{graphs}
\keyword{models}
