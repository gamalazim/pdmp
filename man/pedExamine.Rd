% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pedFuns.R
\name{pedExamine}
\alias{pedExamine}
\title{Examine pedigree data for specific inconsistencies
\code{pedExamine} receives a messy pedigree and the column locations of the individual and
the two parents. The function then performs a specific integrity check and output results as
a list of problematic individuals.}
\usage{
pedExamine(p, cfm = 1:3, examine = c("duplicate", "parent.record",
  "overlap", "misid"), signal = F)
}
\arguments{
\item{p}{data frame or matrix with at least 3 columns for the individual and two parents.}

\item{cfm}{vector of column names or positions of individual, paternal parent and maternal parent.}

\item{recode}{logical for whether or not to output the sorted pedigree in integer-coded format.}
}
\value{
data frame for the sorted pedigree
}
\description{
Examine pedigree data for specific inconsistencies
\code{pedExamine} receives a messy pedigree and the column locations of the individual and
the two parents. The function then performs a specific integrity check and output results as
a list of problematic individuals.
}
\details{
Given an input of consistent pedigree, chronological odering is performed on the pedigree.
The input pedigree is checked first by \code{pedExamine} for all inconsistencie and an error is
generated if inconsistencies are encountered. It is recommended to always run \code{pedCensor} before
\code{pedSort}. See \sQuote{Examples} below.
}
\examples{
## Load pedigree data
data(pedigree)

## Restructure pedigree data in 'p01' and remove inconsistencies
p.clean <- pedCensor(p01)

## Sort in chronological order
p.sort <- pedSort(p.clean)

## Sort pedigree and generate integer-coded output
p.sort.recoded <- pedSort(p.org = p.clean, recode = TRUE)

}
\seealso{
\code{\link{pedCensor}}, \code{\link{pedExamine}}
}

