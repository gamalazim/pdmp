% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pedFuns.R
\name{pedRecode}
\alias{pedRecode}
\title{Recode pedigree data using integer numbers}
\usage{
pedRecode(p, cfm = 1:3)
}
\arguments{
\item{p}{data frame or matrix with at least 3 columns for the individual and two parents.}

\item{cfm}{vector of column names or positions of individual, paternal parent and maternal parent.}
}
\value{
data frame for the recoded pedigree where columns \code{p[,cfm]} are replaced by integer
values consistent across the 3 columns.
}
\description{
This function expects a pedigree and the column locations of the individual and the two parents.
The function then transforms the character IDs in the 3 columns into integer values. The function
recods IDs consistantly across the 3 columns so that a parent that also appears as an individual
will have a single integer ID.
}
\details{
The column of individuals determined by \code{p[,cfm[1]]} is coded from 1 to N. Any parent
who appears as an individual will take his/her integer value of the individuals' column. This
process of recoding is normally performed after structuring and sorting the pedigree. Usually
pedigree data would be free of any inconsistencies before performing this recoding. See 'Examples'.
}
\examples{
## Load pedigree data
data(pedigree)

## Restructure pedigree data in 'p01' and remove inconsistencies
p.clean <- pedCensor(p01)

## Sort in chronological order
p.sort <- pedSort(p.clean)

## Recode sorted pedigree
p.recode <- pedRecode(p.sort)

}
\seealso{
\code{\link{pedCensor}}, \code{\link{pedSort}}
}

