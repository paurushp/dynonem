\name{sim.intervention}
\alias{sim.intervention}
\alias{sim.intervention}

\title{Simulate interventions in a Nested Effects Model}
\description{   
 Simulates a knock-down of a list of network nodes and returns the network nodes and effect reporters, where effects are expected.
}
\usage{
sim.intervention(x, int)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{nem object}
  \item{int}{a character vector of nodes in the network}  
}
\value{
list with two slots:
\item{Sgenes.effected}{effected network nodes}
\item{Egenes.effected}{effected downstream effect reporters}
}

\author{Holger Froehlich}

\examples{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]
   res = nem(D, control=set.default.parameters(unique(colnames(D)), para=c(0.13,0.05))) 
   sim.intervention(res, "rel") # simulate knock-down of rel
}

\keyword{models}
