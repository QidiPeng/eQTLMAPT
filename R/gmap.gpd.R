#'
#' @title Genomic Mediation analysis with Fixed Permutation scheme and
#'   Generalized Pareto Distribution(GPD)
#'
#' @description The gmap.gpd function performs genomic mediation analysis with
#'   Fixed Permutation scheme. It tests for mediation effects for a set of user
#'   specified mediation trios(e.g., eQTL, cis- and trans-genes) in the genome
#'   with the assumption of the presence of cis-association. When the empirical
#'   P-value is small enough, the GPD fit is used to estimate a more accurate
#'   empirical P value.
#'
#'   It returns the mediation p-values(nominal P-value, empirical P-values
#'   obtained from ordinary calculations and empirical P-values estimated using
#'   GPD fitting), the coefficient of linear models(e.g, t_stat, std.error,
#'   beta, beta.total) and the proportions mediated(e.g., the percentage of
#'   reduction in trans-effects after accounting for cis-mediation).
#'
#' @details The function performs genomic mediation analysis with Fixed
#'   Permutation scheme. \code{Adaptive Permutation scheme}{When using Fixed
#'   Permutation scheme, good estimation of insignificant adjusted P-values can
#'   be achieved with few permutations while many more are needed to estimate
#'   highly significant ones. Therefore, we implemented an alternative
#'   permutation scheme that adapts the number of permutations to the
#'   significance level of the variant–phenotype pairs} \code{calculate
#'   Empirical P-values using GPD fitting}{The use of a fixed number of
#'   permutations to calculate empirical P-values has the disadvantage that the
#'   minimum empirical P-value that can be calculated is 1/N. This makes a
#'   larger number of permutations needed to calculate a smaller P-value.
#'   Therefore, we model the tail of the permutation value as a Generalized
#'   Pareto Distribution(GPD), enabling a smaller empirical P-value with fewer
#'   permutation times.}
#'
#' @param snp.dat The eQTL genotype matrix. Each row is an eQTL, each column is
#'   a sample.
#' @param fea.dat A feature profile matrix. Each row is for one feature, each
#'   column is a sample.
#' @param conf A confounders matrix which is adjusted in all mediation tests.
#'   Each row is a confounder, each column is a sample.
#' @param trios.idx A matrix of selected trios indexes (row numbers) for
#'   mediation tests. Each row consists of the index (i.e., row number) of the
#'   eQTL in \code{snp.dat}, the index of cis-gene feature in \code{fea.dat},
#'   and the index of trans-gene feature in \code{fea.dat}. The dimension is the
#'   number of trios by three.
#' @param cl Parallel backend if it is set up. It is used for parallel
#'   computing. We set \code{cl}=NULL as default.
#' @param Minperm The minimum number of permutations. When the number of
#'   permutation statistics better than the original statistic is greater than
#'   \code{Minperm}, stop permutation and directly calculate the empirical P
#'   value. If \code{Minperm}=0, only the nominal P-value is calculated. We set
#'   \code{Minperm}=100 as default.
#' @param Maxperm Maximum number of permutation. We set \code{Maxperm}=10000 as
#'   default.
#' @param gpd.perm Decide when to use GPD to fit estimation parameters. When the
#'   proportion of permutation better than the original statistic is greater
#'   than par, the GPD is fitted to estimate the empirical P-value. We set
#'   \code{gpd.perm}=0.01 as default.
#'
#' @return The algorithm will return a list of nperm, empirical.p,
#'   empirical.p.gpd, nominal.p, beta, std.error, t_stat, beta.total,
#'   beta.change. \item{nperm}{The actual number of permutations for testing
#'   mediation.} \item{empirical.p}{The mediation empirical P-values with nperm
#'   times permutation. A matrix with dimension of the number of trios.}
#'   \item{empirical.p.gpd}{The mediation Empirical P-values with nperm times
#'   permutation using GPD fit. A matrix with dimension of the number of trios.}
#'   \item{nominal.p}{The mediation nominal P-values. A matrix with dimension of
#'   the number of trios.} \item{std.error}{The return std.error value of
#'   feature1 for fit liner models. A matrix with dimension of the number of
#'   trios.} \item{t_stat}{The return t_stat value of feature1 for fit liner
#'   models. A matrix with dimension of the number of trios.} \item{beta}{The
#'   return beta value of feature2 for fit liner models in the case of feature1.
#'   A matrix with dimension of the number of trios.} \item{beta.total}{The
#'   return beta value of feature2 for fit liner models without considering
#'   feature1. A matrix with dimension of the number of trios.}
#'   \item{beta.change}{The proportions mediated. A matrix with dimension of the
#'   number of trios.}
#'
#' @references Ongen H, Buil A, Brown AA, Dermitzakis ET, Delaneau O. (2016)
#'   Fast and efficient QTL mapper for thousands of molecular phenotypes.
#'   Bioinformatics. 2016;32:1479–1485. \doi{10.1093/bioinformatics/btv722}
#' @references Knijnenburg TA, Wessels LFA, Reinders MJT, Shmulevich I. (2009)
#'   Fewer permutations, more accurate P-values. Bioinformatics.
#'   2009;25:i161–i168. \doi{10.1093/bioinformatics/btp211}
#'
#' @examples
#'
#' output <- gmap.gpd(conf = dat$known.conf, fea.dat = dat$fea.dat, snp.dat = dat$snp.dat,
#'                trios.idx = dat$trios.idx[1:10,], Minperm = 10, Maxperm = 1000)
#'
#' \dontrun{
#'   ## generate a cluster with 2 nodes for parallel computing
#'   cl <- makeCluster(2)
#'
#'   ## When the empirical P-value is less than 0.02, a more accurate
#'      empirical P-value is estimated using the GPD fit.
#'   output <- gmap.gpd(conf = dat$known.conf, fea.dat = dat$fea.dat,
#'                      snp.dat = dat$snp.dat, trios.idx = dat$trios.idx[1:10,],
#'                      cl = cl, Minperm = 10, Maxperm = 1000, gpd.perm = 0.02)
#'
#'   stopCluster(cl)
#' }
#'
#' @export
#' @importFrom parallel parLapply
#'
gmap.gpd <- function(snp.dat, fea.dat, conf, trios.idx, cl = NULL, Minperm = 100, Maxperm = 10000, gpd.perm = 0.01){
	confounders <- t(conf)

	triomatrix <- array(NA, c(dim(fea.dat)[2], dim(trios.idx)[1], 3))
	for (i in 1:dim(trios.idx)[1]) {
		triomatrix[,i, ] <- cbind(round(snp.dat[trios.idx[i, 1], ], digits = 0),
								  fea.dat[trios.idx[i, 2], ], fea.dat[trios.idx[i, 3], ])
	}

	num_trio <- dim(triomatrix)[2]

	if(!is.null(cl)){
		output <- parLapply(cl, 1:num_trio, getp.func, triomatrix = triomatrix, confounders = confounders,
		                    Minperm = Minperm, Maxperm = Maxperm, use.gpd=TRUE, gpd.perm = gpd.perm)
	}else{
		output <- lapply(1:num_trio, getp.func, triomatrix = triomatrix, confounders = confounders,
		                 Minperm = Minperm, Maxperm = Maxperm, use.gpd=TRUE, gpd.perm = gpd.perm)
	}

	nominal.p <- matrix(unlist(lapply(output, function(x) x$nominal.p), use.names = FALSE), byrow = T, ncol = 1)
	t_stat <- matrix(unlist(lapply(output, function(x) x$t_stat), use.names = FALSE), byrow = T, ncol = 1)
	std.error <- matrix(unlist(lapply(output, function(x) x$std.error), use.names = FALSE), byrow = T, ncol = 1)
	beta <- matrix(unlist(lapply(output, function(x) x$beta), use.names = FALSE), byrow = T, ncol = 1)
	beta.total <- matrix(unlist(lapply(output, function(x) x$beta.total), use.names = FALSE), byrow = T, ncol = 1)
	beta.change <- matrix(unlist(lapply(output, function(x) x$beta.change), use.names = FALSE), byrow = T, ncol = 1)
	empirical.p <- matrix(unlist(lapply(output, function(x) x$empirical.p), use.names = FALSE), byrow = T, ncol = 1)
	empirical.p.gpd <- matrix(unlist(lapply(output, function(x) x$empirical.p.gpd), use.names = FALSE), byrow = T, ncol = 1)
	nperm <- matrix(unlist(lapply(output, function(x) x$nperm), use.names = FALSE), byrow = T, ncol = 1)
	runtime <- matrix(unlist(lapply(output, function(x) x$runtime), use.names = FALSE), byrow = T, ncol = 1)

	output <- list(nperm = nperm, empirical.p = empirical.p, empirical.p.gpd = empirical.p.gpd, nominal.p = nominal.p,
	               std.error = std.error, t_stat = t_stat, beta = beta, beta.total = beta.total, beta.change = beta.change, runtime = runtime)

	return(output)
}
