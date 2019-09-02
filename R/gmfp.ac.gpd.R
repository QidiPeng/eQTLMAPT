#'
#' @title Genomic Mediation analysis with Fixed Permutation scheme and Adaptive
#'   Confunders and Generalized Pareto Distribution(GPD)
#'
#' @description The gmfp.ac.gpd function performs genomic mediation analysis
#'   with Fixed Permutation scheme and Adaptive Confunders. It tests for
#'   mediation effects for a set of user specified mediation trios(e.g., eQTL,
#'   cis- and trans-genes) in the genome with the assumption of the presence of
#'   cis-association. The gmfp.ac.gpd function considers either a user provided
#'   pool of potential confounding variables, real or constructed by other
#'   methods, or all the PCs based on the feature data as the potential
#'   confounder pool. When the empirical P-value is small enough, the GPD fit is
#'   used to estimate a more accurate empirical P value.
#'
#'   It returns the mediation p-values(nominal P-value, empirical P-values
#'   obtained from ordinary calculations and empirical P-values estimated using
#'   GPD fitting), the coefficient of linear models(e.g, t_stat, std.error,
#'   beta, beta.total), and the proportions mediated(e.g., the percentage of
#'   reduction in trans-effects after accounting for cis-mediation) based on the
#'   mediation tests i) adjusting for known confounders only, and ii) adjusting
#'   for known confounders and adaptively selected potential confounders for
#'   each mediation trio.
#'
#' @details The function performs genomic mediation analysis with Fixed
#'   Permutation scheme and Adaptive Confunders. \code{Fixed Permutation
#'   scheme}{When calculating the empirical P-value, the data is permutated by a
#'   fixed number of times, and the statistics after permutation are separately
#'   calculated. Assuming that the number of permutation is N, where the number
#'   of permutation statistics that is better than the original statistic is M,
#'   then the Empirical P-value = (M + 1) / (N + 1).} \code{Adaptive Confunding
#'   adjustment} {One challenge in mediation test in genomic studies is how to
#'   adjust unmeasured confounding variables for the cis- and trans-genes (i.e.,
#'   mediator-outcome) relationship.The current function adaptively selects the
#'   variables to adjust for each mediation trio given a large pool of
#'   constructed or real potential confounding variables. The function allows
#'   the input of variables known to be potential cis- and trans-genes
#'   (mediator-outcome) confounders in all mediation tests (\code{known.conf}),
#'   and the input of the pool of candidate confounders from which potential
#'   confounders for each mediation test will be adaptively selected
#'   (\code{cov.pool}). When no pool is provided (\code{cov.pool = NULL}), all
#'   the PCs based on feature profile (\code{fea.dat}) will be constructed as
#'   the potential confounder pool.} \code{calculate Empirical P-values using
#'   GPD fitting}{The use of a fixed number of permutations to calculate
#'   empirical P-values has the disadvantage that the minimum empirical P-value
#'   that can be calculated is 1/N. This makes a larger number of permutations
#'   needed to calculate a smaller P-value. Therefore, we model the tail of the
#'   permutation value as a Generalized Pareto Distribution(GPD), enabling a
#'   smaller empirical P-value with fewer permutation times.}
#'
#' @param snp.dat The eQTL genotype matrix. Each row is an eQTL, each column is
#'   a sample.
#' @param fea.dat A feature profile matrix. Each row is for one feature, each
#'   column is a sample.
#' @param known.conf A confounders matrix which is adjusted in all mediation
#'   tests. Each row is a confounder, each column is a sample.
#' @param trios.idx A matrix of selected trios indexes (row numbers) for
#'   mediation tests. Each row consists of the index (i.e., row number) of the
#'   eQTL in \code{snp.dat}, the index of cis-gene feature in \code{fea.dat},
#'   and the index of trans-gene feature in \code{fea.dat}. The dimension is the
#'   number of trios by three.
#' @param cl Parallel backend if it is set up. It is used for parallel
#'   computing. We set \code{cl}=NULL as default.
#' @param cov.pool The pool of candidate confounding variables from which
#'   potential confounders are adaptively selected to adjust for each trio. Each
#'   row is a covariate, each column is a sample. We set \code{cov.pool}=NULL as
#'   default, which will calculate PCs of features as cov.pool.
#' @param pc.num If \code{cov.pool}=NULL, use the previous num PCs as 
#'   \code{cov.pool}.We set \code{pc.num}=30 as default.
#' @param nperm The number of permutations for testing mediation. If
#'   \code{nperm}=0, only the nominal P-value is calculated. We set
#'   \code{nperm}=10000 as default.
#' @param gpd.perm Decide when to use GPD to fit estimation parameters. When the
#'   proportion of permutation better than the original statistic is greater
#'   than par, the GPD is fitted to estimate the empirical P-value. We set
#'   \code{gpd.perm}=0.01 as default.
#' @param fdr The false discovery rate to select confounders. We set
#'   \code{fdr}=0.05 as default.
#' @param fdr_filter The false discovery rate to filter common child and
#'   intermediate variables. We set \code{fdr_filter}=0.1 as default.
#'
#' @return The algorithm will return a list of empirical.p, empirical.p.gpd,
#'   nominal.p, beta, std.error, t_stat, beta.total, beta.change.
#'   \item{empirical.p}{The mediation Empirical P-values with nperm times
#'   permutation. A matrix with dimension of the number of trios.}
#'   \item{empirical.p.gpd}{The mediation empirical P-values with nperm times
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
#'   number of trios.} \item{pc.matrix}{PCs will be returned if the PCs based on
#'   feature data are used as the pool of potential confounders. Each column is
#'   a PC.} \item{sel.conf.ind}{An indicator matrix with dimension of the number
#'   of trios by the number of covariates in \code{cov.pool} or
#'   \code{pc.matrix}if the principal components (PCs) based on feature data are
#'   used as the pool of potential confounders.}
#'
#' @references Yang F, Wang J, Consortium G, Pierce BL, Chen LS. (2017)
#'   Identifying cis-mediators for trans-eQTLs across many human tissues using
#'   genomic mediation analysis. Genome Research. 2017;27:1859–1871.
#'   \doi{10.1101/gr.216754.116}
#' @references Knijnenburg TA, Wessels LFA, Reinders MJT, Shmulevich I. (2009)
#'   Fewer permutations, more accurate P-values. Bioinformatics.
#'   2009;25:i161–i168. \doi{10.1093/bioinformatics/btp211}
#'
#' @examples
#'
#' output <- gmfp.ac.gpd(known.conf = dat$known.conf, fea.dat = dat$fea.dat,
#'                       snp.dat = dat$snp.dat, trios.idx = dat$trios.idx[1:10,], nperm = 100)
#'
#' \dontrun{
#'   ## generate a cluster with 2 nodes for parallel computing
#'   cl <- makeCluster(2)
#'
#'   ## Use the specified candidate confusion variable pool
#'   ## When the empirical P-value is less than 0.02, a more accurate
#'      empirical P-value is estimated using the GPD fit.
#'   output <- gmfp.ac.gpd(known.conf = dat$known.conf, fea.dat = dat$fea.dat,
#'                         snp.dat = dat$snp.dat, trios.idx = dat$trios.idx[1:10,],
#'                         cl = cl, cov.pool = dat$cov.pool, nperm = 100, gpd.perm = 0.02)
#'
#'   stopCluster(cl)
#' }
#'
#' @export
#' @importFrom parallel parLapply
#'
gmfp.ac.gpd <- function(snp.dat, fea.dat, known.conf, trios.idx, cl = NULL, cov.pool = NULL, pc.num = 30,
						nperm = 10000, gpd.perm = 0.01, fdr = 0.05, fdr_filter = 0.1){
	confounders <- t(known.conf)

	triomatrix <- array(NA, c(dim(fea.dat)[2], dim(trios.idx)[1], 3))
	for (i in 1:dim(trios.idx)[1]) {
		triomatrix[,i, ] <- cbind(round(snp.dat[trios.idx[i, 1], ], digits = 0),
								  fea.dat[trios.idx[i, 2], ], fea.dat[trios.idx[i, 3], ])
	}

	num_trio <- dim(triomatrix)[2]

	# Adaptive Confunding adjustment
	use.PC <- FALSE
	if(is.null(cov.pool)){ # using PC
		use.PC <- TRUE
		res <- get.cov(cl, fea.dat = fea.dat, triomatrix = triomatrix, fdr = fdr, fdr_filter = fdr_filter)
		pool_cov <- res$pool_cov
		est_conf_pool_idx <- res$est_conf_pool_idx
		all_pc <- res$pc.all
	}else{ # using cov.pool
		res <- get.cov(cl, cov.pool = cov.pool, triomatrix = triomatrix, fdr = fdr, fdr_filter = fdr_filter)
		pool_cov <- res$pool_cov
		est_conf_pool_idx <- res$est_conf_pool_idx
	}
	est_conf_pool_idx <- est_conf_pool_idx[ ,c(1:pc.num)]

	if(!is.null(cl)){
		known_output <- parLapply(cl, 1:num_trio, getp.func, triomatrix = triomatrix, confounders = confounders,
		                          Minperm = nperm, Maxperm = nperm, use.gpd = TRUE, gpd.perm = gpd.perm)
		known_sel_pool_output <- parLapply(cl, 1:num_trio, getp.func, triomatrix = triomatrix,
										   confounders = confounders, pool_cov = pool_cov, est_conf_pool_idx = est_conf_pool_idx,
										   Minperm = nperm, Maxperm = nperm, use.gpd = TRUE, gpd.perm = gpd.perm, use.PC = use.PC)
	}else{
		known_output <- lapply(1:num_trio, getp.func, triomatrix = triomatrix, confounders = confounders,
		                       Minperm = nperm, Maxperm = nperm, use.gpd = TRUE, gpd.perm = gpd.perm)
		known_sel_pool_output <- lapply(1:num_trio, getp.func, triomatrix = triomatrix,
						 				confounders = confounders, pool_cov = pool_cov, est_conf_pool_idx = est_conf_pool_idx,
						 				Minperm = nperm, Maxperm = nperm, use.gpd = TRUE, gpd.perm = gpd.perm, use.PC = use.PC)
	}

	nominal.p <- matrix(c(lapply(known_output, function(x) x$nominal.p),
	                      lapply(known_sel_pool_output, function(x) x$nominal.p)), byrow = F, ncol = 2)
	t_stat <- matrix(c(lapply(known_output, function(x) x$t_stat),
	                   lapply(known_sel_pool_output, function(x) x$t_stat)), byrow = F, ncol = 2)
	std.error <- matrix(c(lapply(known_output, function(x) x$std.error),
	                      lapply(known_sel_pool_output, function(x) x$std.error)), byrow = F, ncol = 2)
	beta <- matrix(c(lapply(known_output, function(x) x$beta),
	                 lapply(known_sel_pool_output, function(x) x$beta)), byrow = F, ncol = 2)
	beta.total <- matrix(c(lapply(known_output, function(x) x$beta.total),
							lapply(known_sel_pool_output, function(x) x$beta.total)), byrow = F, ncol = 2)
	beta.change <- matrix(c(lapply(known_output, function(x) x$beta.change),
							lapply(known_sel_pool_output, function(x) x$beta.change)), byrow = F, ncol = 2)
	empirical.p <- matrix(c(lapply(known_output, function(x) x$empirical.p),
							lapply(known_sel_pool_output, function(x) x$empirical.p)), byrow = F, ncol = 2)
	empirical.p.gpd <- matrix(c(lapply(known_output, function(x) x$empirical.p.gpd),
								lapply(known_sel_pool_output, function(x) x$empirical.p.gpd)), byrow = F, ncol = 2)
#	nperm <- matrix(c(lapply(known_output, function(x) x$nperm),
#	                  lapply(known_sel_pool_output, function(x) x$nperm)), byrow = F, ncol = 2)

	if(use.PC){
		output <- list(empirical.p = empirical.p, empirical.p.gpd = empirical.p.gpd, nominal.p = nominal.p,
		               std.error = std.error, t_stat = t_stat, beta = beta, beta.total = beta.total, beta.change = beta.change,
		               pc.matrix = all_pc, sel.conf.ind = est_conf_pool_idx)
	}else{
		output <- list(empirical.p = empirical.p, empirical.p.gpd = empirical.p.gpd, nominal.p = nominal.p,
		               std.error = std.error, t_stat = t_stat, beta = beta, beta.total = beta.total, beta.change = beta.change,
		               sel.conf.ind = est_conf_pool_idx)
	}

	return(output)
}
