#'
#' @title Genomic Mediation analysis with Fixed Permutation scheme
#'
#' @description The gmfp function performs genomic mediation analysis with fixed
#'   permutation. It tests for mediation effects for a set of user specified
#'   mediation trios(e.g., eQTL, cis- and trans-genes) in the genome with the
#'   assumption of the presence of cis-association.
#'
#'   It returns the mediation p-values(nominal and empirical), the coefficient
#'   of linear models(e.g, t_stat, std.error, beta, beta.total) and the
#'   proportions mediated(e.g., the percentage of reduction in trans-effects
#'   after accounting for cis-mediation).
#'
#' @details The function performs genomic mediation analysis with fixed
#'   permutation. \code{Fixed Permutation scheme}{When calculating the empirical
#'   P-value, the data is permutated by a fixed number of times, and the
#'   statistics after permutation are separately calculated. Assuming that the
#'   number of permutation is N, where the number of permutation statistics that
#'   is better than the original statistic is M, then the Empirical P-value = (M
#'   + 1) / (N + 1).}
#'
#' @param snp.dat The eQTL genotype matrix. Each row is an eQTL, each column is
#'   a sample.
#' @param fea.dat A feature profile matrix. Each row is for one feature, each
#'   column is a sample.
#' @param conf A confounders matrix which is adjusted in mediation tests. Each
#'   row is a confounder, each column is a sample.
#' @param trios.idx A matrix of selected trios indexes (row numbers) for
#'   mediation tests. Each row consists of the index (i.e., row number) of the
#'   eQTL in \code{snp.dat}, the index of cis-gene feature in \code{fea.dat},
#'   and the index of trans-gene feature in \code{fea.dat}. The dimension is the
#'   number of trios by three.
#' @param cl Parallel backend if it is set up. It is used for parallel
#'   computing. We set \code{cl}=NULL as default.
#' @param nperm The number of permutations for testing mediation. If
#'   \code{nperm}=0, only the nominal P-value is calculated. We set
#'   \code{nperm}=10000 as default.
#'
#' @return The algorithm will return a list of empirical.p, nominal.p, beta,
#'   std.error, t_stat, beta.total, beta.change. \item{empirical.p}{The
#'   mediation empirical P-values with nperm times permutation. A matrix with
#'   dimension of the number of trios.} \item{nominal.p}{The mediation nominal
#'   P-values. A matrix with dimension of the number of trios.}
#'   \item{std.error}{The return std.error value of feature1 for fit liner
#'   models. A matrix with dimension of the number of trios.} \item{t_stat}{The
#'   return t_stat value of feature1 for fit liner models. A matrix with
#'   dimension of the number of trios.} \item{beta}{The return beta value of
#'   feature2 for fit liner models in the case of feature1. A matrix with
#'   dimension of the number of trios.} \item{beta.total}{The return beta value
#'   of feature2 for fit liner models without considering feature1. A matrix
#'   with dimension of the number of trios.} \item{beta.change}{The proportions
#'   mediated. A matrix with dimension of the number of trios.}
#'
#' @examples
#'
#' output <- gmfp(conf = dat$known.conf, fea.dat = dat$fea.dat, snp.dat = dat$snp.dat,
#'                trios.idx = dat$trios.idx[1:10,], nperm = 100)
#'
#' \dontrun{
#'   ## generate a cluster with 2 nodes for parallel computing
#'   cl <- makeCluster(2)
#'   output <- gmfp(conf = dat$known.conf, fea.dat = dat$fea.dat, snp.dat = dat$snp.dat,
#'                  trios.idx = dat$trios.idx[1:10,], cl = cl, nperm = 100)
#'   stopCluster(cl)
#' }
#'
#' @export
#' @importFrom parallel parLapply
#'
gmfp <- function(snp.dat, fea.dat, conf, trios.idx, cl = NULL, nperm = 10000){
	confounders <- t(conf)

	triomatrix <- array(NA, c(dim(fea.dat)[2], dim(trios.idx)[1], 3))
	for (i in 1:dim(trios.idx)[1]) {
	    triomatrix[,i, ] <- cbind(round(snp.dat[trios.idx[i, 1], ], digits = 0),
	                              fea.dat[trios.idx[i, 2], ], fea.dat[trios.idx[i, 3], ])
	}

	num_trio <- dim(triomatrix)[2]

	if(!is.null(cl)){
		output <- parLapply(cl, 1:num_trio, getp.func, triomatrix = triomatrix, confounders = confounders,
								Minperm = nperm, Maxperm = nperm)
	}else{
		output <- lapply(1:num_trio, getp.func, triomatrix = triomatrix, confounders = confounders,
							Minperm = nperm, Maxperm = nperm)
	}

	nominal.p <- matrix(unlist(lapply(output, function(x) x$nominal.p), use.names = FALSE), byrow = T, ncol = 1)
	t_stat <- matrix(unlist(lapply(output, function(x) x$t_stat), use.names = FALSE), byrow = T, ncol = 1)
	std.error <- matrix(unlist(lapply(output, function(x) x$std.error), use.names = FALSE), byrow = T, ncol = 1)
	beta <- matrix(unlist(lapply(output, function(x) x$beta), use.names = FALSE), byrow = T, ncol = 1)
	beta.total <- matrix(unlist(lapply(output, function(x) x$beta.total), use.names = FALSE), byrow = T, ncol = 1)
	beta.change <- matrix(unlist(lapply(output, function(x) x$beta.change), use.names = FALSE), byrow = T, ncol = 1)
	empirical.p <- matrix(unlist(lapply(output, function(x) x$empirical.p), use.names = FALSE), byrow = T, ncol = 1)
#	nperm <- matrix(unlist(lapply(output, function(x) x$nperm), use.names = FALSE), byrow = T, ncol = 1)

	output <- list(empirical.p = empirical.p, nominal.p = nominal.p, std.error = std.error,
	               t_stat = t_stat, beta = beta, beta.total = beta.total, beta.change = beta.change)

	return(output)
}
