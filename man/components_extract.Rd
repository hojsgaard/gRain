% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpt_pot_extract.R
\name{components_extract}
\alias{components_extract}
\alias{extract_cpt}
\alias{extract_pot}
\alias{extract_marg}
\alias{marg2pot}
\alias{pot2marg}
\title{Extract conditional probabilities and clique potentials from
data.}
\usage{
extract_cpt(data_, graph, smooth = 0)

extract_pot(data_, graph, smooth = 0)

extract_marg(data_, graph, smooth = 0)

marg2pot(marg_rep)

pot2marg(pot_rep)
}
\arguments{
\item{data_}{A named array or a dataframe.}

\item{graph}{An \code{igraph} object or a list or formula which can be
turned into a \code{igraph} object by calling \code{ug} or
\code{dag}. For \code{extract_cpt}, graph must be/define a DAG while for
\code{extract_pot}, graph must be/define undirected triangulated graph.}

\item{smooth}{See 'details' below.}

\item{marg_rep}{An object of class \code{marg_rep}}

\item{pot_rep}{An object of class \code{pot_representation}}
}
\value{
\itemize{
\item \code{extract_cpt}: A list of conditional probability tables.
\item \code{extract_pot}: A list of clique potentials.
\item \code{extract_marg}: A list of clique marginals.
}
}
\description{
Extract list of conditional probability tables and
list of clique potentials from data.
}
\details{
If \code{smooth} is non-zero then \code{smooth} is added
to all cell counts before normalization takes place.
}
\examples{

## Extract cpts / clique potentials from data and graph
# specification and create network. There are different ways:

data(lizard, package="gRbase")

# DAG: height <- species -> diam
daG <- dag(~species + height:species + diam:species, result="igraph")

# UG : [height:species][diam:species]
uG  <- ug(~height:species + diam:species, result="igraph")

pt <- extract_pot(lizard, ~height:species + diam:species) 
cp <- extract_cpt(lizard, ~species + height:species + diam:species)

pt
cp

# Both specify the same probability distribution
tabListMult(pt) |> as.data.frame.table()
tabListMult(cp) |> as.data.frame.table()

\dontrun{
# Bayesian networks can be created as
bn.uG   <- grain(pt)
bn.daG  <- grain(cp)

# The steps above are wrapped into a convenience method which
# builds a network from at graph and data.
bn.uG   <- grain(uG, data=lizard)
bn.daG  <- grain(daG, data=lizard)
}

}
\references{
Søren Højsgaard (2012). Graphical Independence Networks
with the gRain Package for R. Journal of Statistical Software,
46(10), 1-26.  \url{https://www.jstatsoft.org/v46/i10/}.
}
\seealso{
\code{\link{compileCPT}}, \code{\link{compilePOT}},
\code{\link{grain}}
}
\author{
Søren Højsgaard, \email{sorenh@math.aau.dk}
}
\keyword{utilities}
