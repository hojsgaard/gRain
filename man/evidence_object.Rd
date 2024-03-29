% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evidence_objects.R
\name{evidence_object}
\alias{evidence_object}
\alias{new_evi}
\alias{subset.grain_evidence}
\alias{print.grain_evidence}
\alias{varNames.grain_evidence}
\alias{is.null_evi}
\alias{as.data.frame.grain_evidence}
\alias{setdiff_evi}
\alias{union_evi}
\title{Evidence objects}
\usage{
new_evi(evi_list = NULL, levels)

is.null_evi(object)

\method{as.data.frame}{grain_evidence}(x, row.names = NULL, optional = FALSE, ...)

setdiff_evi(ev1, ev2)

union_evi(ev1, ev2)
}
\arguments{
\item{evi_list}{A named list with evidence; see 'examples' below.}

\item{levels}{A named list with the levels of all variables.}

\item{object}{Some R object.}

\item{x}{An evidence object.}

\item{row.names}{Not used.}

\item{optional}{Not used.}

\item{...}{Not used.}

\item{ev1, ev2}{Evidence.}
}
\description{
Functions for defining and manipulating evidence.
}
\details{
Evidence is specified as a list. Internally, evidence is
represented as a grain evidence object which is a list with 4 elements.
}
\examples{

## Define the universe
yn <- c("yes", "no")
uni <- list(asia = yn, tub = yn, smoke = yn, lung = yn,
            bronc = yn, either = yn, xray = yn, dysp = yn)

e1 <- list(dysp="no", xray="no")
eo1 <- new_evi(e1, levels=uni)
eo1  |> as.data.frame()

e2 <- list(dysp="no", xray=c(0, 1))
eo2 <- new_evi(e2, levels=uni)
eo2 |> as.data.frame()

# Above e1 and e2 specifies the same evidence but information about
# whether the state has been set definite or as a weight is
# maintained.

e3 <- list(dysp="yes", asia="yes")
eo3 <- new_evi(e3, uni)
eo3 |> as.data.frame()

# If evidence 'e1' is already set in the network and new evidence
# 'e3' emerges, then evidence in the network must be updated. But
# there is a conflict in that dysp="yes" in 'e1' and
# dysp="no" in 'e3'. The (arbitrary) convention is that
# existing evidence overrides new evidence so that the only new
# evidence in 'e3' is really asia="yes".

# To subtract existing evidence from new evidence we can do:
setdiff_evi(eo3, eo1) |> as.data.frame()

# Likewise the 'union' is
union_evi(eo3, eo1) |> as.data.frame()

}
\author{
Søren Højsgaard, \email{sorenh@math.aau.dk}
}
