
get.lm <- function(mediator, treatment, outcome, confounderset, return.t.only = TRUE){
	z <- summary(lm(outcome ~ mediator + treatment + confounderset))
	t_stat <- z$coef[,3][2]

	if(return.t.only)
		return(t_stat = t_stat)
	else{
		pval <- z$coef[,4][2]
		std <- z$coef[,2][2]
		z2 <- summary(lm(outcome ~ treatment + confounderset))
		beta <- z$coef[,1][3]
		beta.total <- z2$coef[,1][2]
		return(list(pval = pval, t_stat = t_stat, std = std, beta = beta, beta.total = beta.total))
	}
}

get.beta.change <- function(beta_direct, beta_total) {
	beta_change <- (beta_total - beta_direct) / beta_total
	return(beta_change)
}


#' @title This function caculate P-value for every trio
#' @description This function calculates the P-value for each trios. If
#'   \code{Minperm}=0, only the nominal P-value is calculated. If
#'   \code{Minperm}=\code{Maxperm}, the empirical P-value is calculated using a
#'   fixed number of permutation statistics; otherwise, the empirical P-value is
#'   calculated using the adaptive permutation scheme. The user can specify
#'   whether to use the GPD fit to estimate a more accurate empirical P-value,
#'   and at the same time specify how small the empirical P-value is for GPD
#'   fitting.
#'
#' @param i Trios index in triomatrix
#' @param triomatrix A three-dimensional matrix of size: samples number * trios
#'   number * 3. Triomatrix[i,j,1] represents the genotype of the j-th trios at
#'   the i-th sample,and triomatrix[i,j,2] represents the feature1 data of the
#'   j-th trios at the i-th sample, triomatrix[i,j,3] represents the feature2
#'   data of the j-th trios at the i-th sample.
#' @param confunders A confounders matrix which is adjusted in all mediation
#'   tests.
#' @param Minperm Decide whether to use the parameters of the GPD fit. If the
#'   value is 0, only the nominal P-value is calculated. If the proportion of
#'   the permutation statistic better than the original statistic to the total
#'   number of permutations exceeds this value, a more accurate empirical P
#'   value is estimated using the GPD fit. If \code{Minperm}=0, only the nominal
#'   P-value is calculated. We set \code{Minperm}=100 as default.
#' @param Maxperm Maximum number of permutation. We set \code{Maxperm}=10000 as
#'   default.
#' @param use.gpd Whether to use the GPD fit to estimate a more accurate
#'   empirical P-value. We set \code{use.gpd}=NULL as default.
#' @param gpd.perm The proportion parameter for estimating the empirical P-value
#'   when using GPD fit. When the proportion of permutation better than the
#'   original statistic is greater than par, the GPD is fitted to estimate the
#'   empirical P-value. We set \code{gpd.perm}=0.01 as default.
#' @param pool_cov Candidate Confusion Variable Pool. We set
#'   \code{pool_cov}=NULL as default, which \code{use.PC} requied.
#' @param est_conf_pool_idx The index of the adaptively selected confunding
#'   variable. We set \code{pool_cov}=NULL as default, at this time
#'   \code{pool_cov}=NULL too.
#' @param use.PC Whether the candidate confusion variable pool is PCs.
#'
#' @return The algorithm will return a list of nperm, empirical.p,
#'   empirical.p.gpd, nominal.p, std.error, t_stat, beta, beta.total,
#'   beta.change. \item{nperm}{The actual number of permutations for testing
#'   mediation, equal to the input parameter nperm.} \item{empirical.p}{The
#'   mediation Empirical P-values with nperm times permutation.}
#'   \item{empirical.p.gpd}{The mediation Empirical P-values with nperm times
#'   permutation using GPD fit.} \item{nominal.p}{The mediation nominal
#'   P-values. A matrix with dimension of the number of trios.}
#'   \item{std.error}{The return std.error value of feature1 for fit liner
#'   models.} \item{t_stat}{The return t_stat value of feature1 for fit liner
#'   models.} \item{beta}{The return beta value of feature2 for fit liner models
#'   in the case of feature1.} \item{beta.total}{The return beta value of
#'   feature2 for fit liner models without considering feature1.}
#'   \item{beta.change}{The proportions mediated.}
#'
getp.func <- function(i, triomatrix, confounders, Minperm=100, Maxperm=10000, use.gpd = FALSE, gpd.perm = 0.01,
						pool_cov = NULL, est_conf_pool_idx = NULL, use.PC = FALSE) {
    begin.time = proc.time()
    
	treatment <- triomatrix[, i, 1] #L
	mediator <- triomatrix[, i, 2] #C
	outcome <- triomatrix[, i, 3] #T

	# Adaptive Confunding adjustment
	if(!is.null(pool_cov) && !is.null(est_conf_pool_idx)){
		if(use.PC)
			pool_cov <- pool_cov[-nrow(pool_cov),]
		idx = which(est_conf_pool_idx[i, ] == 1)
		if(length(idx) != 0)
	    	confounders <- cbind(confounders, pool_cov[, idx])
	}

	indirect <- get.lm(mediator = mediator, treatment = treatment, outcome = outcome,
						confounderset = confounders, return.t.only = FALSE)
	t_stat <- indirect$t_stat
	nominal.p <- indirect$pval
	beta.change <- get.beta.change(beta_direct = indirect$beta, beta_total = indirect$beta.total)

	## only nominal.p
	if(Minperm == 0){
		output <- list(nominal.p = nominal.p, t_stat = t_stat, std.error = indirect$std,
						beta = indirect$beta, beta.total = indirect$beta.total, beta.change = beta.change, nperm = 0)
		return(output)
	}

	if(Minperm < Maxperm){ # adaptive permutation scheme
	    t_stat_perm <- c()
	    b = 0
	    for(i in 1:Maxperm){
	        mediator_perm <- mediator
	        for(j in 0:2){
	            ind = which(treatment == j)
	            if(length(ind) > 1)
	                mediator_perm[ind] <- sample(mediator_perm[ind], size = length(ind))
	        }
	        
	        t_stat_permutation <- get.lm(mediator_perm, treatment = treatment, outcome = outcome, confounderset = confounders)
	        t_stat_perm <- c(t_stat_perm, t_stat_permutation)
	        if(abs(t_stat) <= abs(t_stat_permutation))
	            b = b + 1
	        if(b >= Minperm)
	            break
	    }
	}else{ # direct permutation scheme
	    mediator_perm <- matrix(rep(mediator, each = Maxperm), byrow = TRUE, ncol = Maxperm)
	    for (j in 0:2) {
	        ind <- which(treatment == j)
	        if (length(ind) > 1){
	            mediator_perm[ind, ] <- apply(mediator_perm[ind, ], 2, sample)
	        }
	    }
	    t_stat_perm <- apply(mediator_perm, 2, get.lm, treatment = treatment, outcome = outcome, confounderset = confounders)
	}

	## empirical.p
	nperm <- length(t_stat_perm)
	r <- sum(abs(t_stat) <= abs(t_stat_perm))
	empirical.p <- (r+1)/(nperm+1)

	end.time = proc.time()
	## empirical.p using GPD fit
	if(use.gpd){
		empirical.p.gpd <- Ppermest(abs(t_stat), abs(t_stat_perm), proportion = gpd.perm)
		output <- list(nperm = nperm, nominal.p = nominal.p, t_stat = t_stat, std.error = indirect$std,
		               beta = indirect$beta, beta.total = indirect$beta.total, beta.change = beta.change,
		               empirical.p = empirical.p, empirical.p.gpd = empirical.p.gpd, runtime = (end.time-begin.time)[1])
	}else{
		output <- list(nperm = nperm, nominal.p = nominal.p, t_stat = t_stat, std.error = indirect$std,
		               beta = indirect$beta, beta.total = indirect$beta.total, beta.change = beta.change, 
		               empirical.p = empirical.p, runtime = (end.time-begin.time)[1])
	}
	return(output)
}

