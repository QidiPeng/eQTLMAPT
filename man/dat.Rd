% Please edit documentation in R/data.R
\docType{data}
\name{dat}
\alias{dat}
\title{Example data}
\value{
A list containing
\item{known.conf}{The known confounders matrix which is adjusted in all mediation tests. Each row is a confounder, each column is a sample.}
\item{cov.pool}{The pool of candidate confounding variables from which potential confounders are adaptively selected to adjust for each mediation test. Each row is a covariate, each column is a sample.}
\item{fea.dat}{The gene expression matrix. Each row is for one gene, each column is a sample.}
\item{snp.dat}{The eQTL genotype matrix. Each row is an eQTL, each column is a sample.}
\item{trios.idx}{The matrix of selected trios indexes (row numbers) for mediation tests. Each row consists of the index (i.e., row number) of the eQTLs in \code{snp.dat},
the index of cis-gene transcript in \code{fea.dat}, and the index of trans-gene transcript in \code{fea.dat}. The dimension is the number of trios by three.}
}
\description{
This simulated data list is for demonstration.
}
\examples{
data(example)
}
