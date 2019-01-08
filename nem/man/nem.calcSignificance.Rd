\name{nem.calcSignificance}
\alias{nem.calcSignificance}

\title{Statistical significance of network hypotheses}
\description{   
 Assess statistical significance of a network hypothesis by comparing it to a null hypothesis. 
}
\usage{
	nem.calcSignificance(D, x, N=1000, seed=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix with experiments in the columns (binary or continious)}  
  \item{x}{nem object}  
  \item{N}{number of random networks to sample}
  \item{seed}{random seed}    
}
\details{
Given data, N random network hypotheses from a null distribution are drawn as follows:  For each S-gene \eqn{S_{k}} we randomly choose a number o of outgoing edges between 0 and 3. We then select o S-genes having at most 1 ingoing edge, connected  \eqn{S_{k}} to them and transitively closed the graph. 
For all random network hypotheses it is counted, how often their likelihood is bigger than that of the given network. This yields an exact p-value.

Another way of assessing the statistical significance of the network hypothesis is to draw random permutations of the node labels. Note that in this case the node degree distribution is the same as in the given network. Again, we can obtain an exact p-value by counting, how often the likelihood of the permuted network is bigger than that of the given network.

Finally, comparison to randomly perturbed networks (insertion or deletion of 1 edge) yields an exact p-value describing the stability of the network.
}
\value{
\item{p.value.rnd}{p-value of the network according to the null hypothesis that it is random}
\item{p.value.perm}{p-value of the network according to the null hypothesis that a network with permuted node labels is at least as good}
\item{p.value.mod}{p-value of the network according to the null hypothesis a randomly peturbed network is at least as good}
}

\author{Holger Froehlich}

\seealso{\code{\link{nem.consensus}}, \code{\link{nem.jackknife}}, \code{\link{nem.bootstrap}}, \code{\link{nem}}}
\examples{
\dontrun{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]   
   res = nem(D, control=set.default.parameters(unique(colnames(D)), para=c(0.13,0.05))) # get best network
   nem.calcSignificance(D,res) # assess its statistical significance
}
}

\keyword{models}
