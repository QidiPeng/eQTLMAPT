## cov.pool 与 fea.dat 不能同时为NULL
#' @title This function get the Adaptive Confunding adjustments
#' @importFrom parallel parLapply
#'
get.cov <- function(cl = NULL, cov.pool = NULL, fea.dat =NULL, triomatrix, fdr = 0.05, fdr_filter = 0.1){
	use.PC <- FALSE
	if(is.null(cov.pool)){
		use.PC <- TRUE
		res <- get.pc(fea.dat)
		cov.pool <- res$cov.pool
		pc.all <- res$all_pc
	}

	pool_cov <- t(cov.pool)
	num_pool <- dim(pool_cov)[2]

	num_trio <- dim(triomatrix)[2]
	## 1.obtain pool confounder candidates and estimate pool confounders
	if(is.null(cl)){
		p_value_child_pool <- matrix(unlist(lapply(1:num_trio, child.p, tripletmatrix = triomatrix, covariates = pool_cov),
											use.names = F), byrow = TRUE, ncol = num_pool)
	}else{
		p_value_child_pool <- matrix(unlist(parLapply(cl, 1:num_trio, child.p,
														tripletmatrix = triomatrix, covariates = pool_cov),
											use.names = F), byrow = TRUE, ncol = num_pool)
	}
	q_child_pool <- qvalue(p = as.vector(p_value_child_pool), fdr.level = fdr_filter)$qvalues
	conf_candidates_pool <- matrix(0, nrow = nrow(p_value_child_pool), ncol = num_pool)
	conf_candidates_pool[which(q_child_pool >= fdr_filter)] <- 1

	## 2.Estimate pool confounders
	if(is.null(cl)){
		est_conf_pool_idx <- matrix(unlist(lapply(1:num_pool, conf.fdr, tripletmatrix = triomatrix,
												  covariates = pool_cov, conf_candidates = conf_candidates_pool, fdr = fdr),
											use.names = F), ncol = num_pool)
	}else{
		est_conf_pool_idx <- matrix(unlist(parLapply(cl, 1:num_pool, conf.fdr, tripletmatrix = triomatrix,
													 covariates = pool_cov, conf_candidates = conf_candidates_pool, fdr = fdr),
											use.names = F), ncol = num_pool)
	}

	if(use.PC){
		output <- list(pc.all = pc.all, est_conf_pool_idx = est_conf_pool_idx, pool_cov = pool_cov)
	}else{
		output <- list(est_conf_pool_idx = est_conf_pool_idx, pool_cov = pool_cov)
	}
	return(output)
}



#' @title normalize fearession data so that every gene contributes equally to
#'   the construction of PCs
get.pc <- function(fea.dat){
	fea.dat2 <- t(apply(fea.dat, 1, function(x) (x - mean(x))/sd(x)))
	pc <- prcomp(t(fea.dat2))
	all_pc <- pc$x
	eigen <- pc$sdev^2
	feavar <- eigen/sum(eigen)
	cov.pool <- t(rbind(all_pc, feavar))
	return(list(cov.pool = cov.pool, all_pc = all_pc))
}


#' @title This function uses stratified fdr to figure out, for each locus, the
#'   list of covariates that do not play roles as child/intermediate mediator.
child.p <- function(i, tripletmatrix, covariates){
	treatment <- tripletmatrix[, i, 1]
	n_obs <- dim(tripletmatrix)[1]
	n_cov <- dim(covariates)[2]
	p_value_child <- rep(NA, n_cov)
	cov_length <- dim(covariates)[1]
	if(cov_length == (n_obs + 1)){
		covariates <- covariates[-cov_length, ]
	}
	for(j in 1:n_cov){
		p_value_child[j] <- summary(lm(covariates[, j] ~ treatment))$coef[2, 4]
	}
	return(p_value_child)
}


#' @title This function uses stratified fdr to figure out, for each covariate,
#'   the list of trios where the covariate plays a role as a confounder.
conf.fdr <- function(i, tripletmatrix, covariates, conf_candidates, fdr){
	n_obs <- dim(tripletmatrix)[1]
	n_tri <- dim(tripletmatrix)[2]
	p_value <- rep(NA, n_tri)
	cov <- covariates[, i]
	cov_length <- length(cov)
	str_fdr <- rep(0, n_tri)
	candidate_trio_id <- sort(which(conf_candidates[, i] == 1))
	if(sum(candidate_trio_id) == 0)
	    return(str_fdr)
	if(length(cov) == (n_obs + 1)){
		pi0 = 1 - cov[cov_length]
		cov <- cov[-cov_length]
		for(j in 1:n_tri){
			f <- summary(lm(cov ~ tripletmatrix[, j, 2] + tripletmatrix[, j, 3]))$fstatistic
			p_value[j] <- pf(f[1], f[2], f[3], lower.tail = F)
		}
		q <- pi0 * p_value[candidate_trio_id] * length(candidate_trio_id)/rank(p_value[candidate_trio_id])
		str_fdr[candidate_trio_id[which(q <= fdr)]] <- 1
	}else if(cov_length == n_obs){
		for(j in 1:n_tri){
			f <- summary(lm(cov ~ tripletmatrix[, j, 2] + tripletmatrix[, j, 3]))$fstatistic
			p_value[j] <- pf(f[1], f[2], f[3], lower.tail = F)
		}
		q <- qvalue(p = p_value[candidate_trio_id], fdr.level = fdr)$qvalues
		str_fdr[candidate_trio_id[which(q <= fdr)]] <- 1
	}
	return(str_fdr)
}

