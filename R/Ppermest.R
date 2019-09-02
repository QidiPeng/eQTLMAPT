#' @title Estimation of P-value using Generalized Pareto Distribution(GPD)
#' @parm x0  test statistic.
#' @parm y  permutation values.
#' @parm Nexcmax  Number of initial exceedances (default: 250).
#' @parm proportion  When to estimate the empirical P value using GPD fitting.
#'   When the number of permutation statistices better than the original
#'   statistic is less than the total number of permutation, the empirical P
#'   value is estimated using the GPD fit.
#' @return P  Estimate of the P-value.
Ppermest <- function(x0, y, Nexcmax = 250, proportion = 0.01){
	y <- sort(y, decreasing=TRUE)
	N <- length(y)
	M <- sum(y >= x0)
	## limiting Nexcmax to a quarter of the distribution, assuming that the tail is smaller than a quarter
	Nexcmax <- min(Nexcmax, N/4)

	if(M >= N * proportion){
		P <- M / N
	}else{
		print("call Pgpd")
		Nexcvec <- seq(Nexcmax, 10, -10)
		LNV <- length(Nexcvec)
		p <- 1
		P <- NaN
		while(is.nan(P) & p <= LNV){
			Nexc <- Nexcvec[p]
			P <- Pgpd(y, x0, N, Nexc)
			p <- p + 1
		}
	}
	return(P)
}



#' @title Computing permutation test P-value of the GPD approximation using 'ML'
#' @parm y  permutation values
#' @parm x0  original statistic
#' @parm N  number of permutation values
#' @parm Nexc  number of permutations used to approximate the tail
#' @return Phat  estimated P-value
Pgpd <- function(y, x0, N, Nexc){
	# Threshold for AD statistic
	Padth <- 0.05

	# Defining the tail
	z <- y[1:Nexc]
	t <- mean(y[Nexc:Nexc+1])
	z <- z - t
	z <- z[z > 0]
	Nexc <- length(z)
	frac <- Nexc / N

	# Fitting the tail and computing the Pvalue
	parmhat <- gpdfit(z) ##make fit
#	a <- parmhat[1]
#	k <- parmhat[2]

	if(!is.nan(parmhat[1]) && !is.nan(parmhat[2])){ ##check if fit was a success
		a <- parmhat[1] # scale parm
		k <- parmhat[2] # shape parm
		res <- gpdgoft(gpdcdf(z,a,k),k) ##goodness-of-fit test
#		Pcm <- res$Pcm
		Pad <- res$Pad
#		W2 <- res$W2
#		A2 <- res$A2

		if(Pad > Padth){ ##check goodness of fit
			## Pgpd = (Nexc/N) * (1-F(x0-t))
			Phat <- (Nexc/N) * gpdcdf(x0-t, parmhat[1], parmhat[2])
		}else
			Phat = NaN
	}else{
		Phat = NaN
	}
	return(Phat)
}



#' @title (1-CDF) of the generalized pareto distribution
#' @parm x  exceedances
#' @parm a  shape parameter
#' @parm k  scale parameter
#' @return p  probability
gpdcdf <- function(z, a, k){
	if(abs(k) == 0) #k==0(<eps)
		p <- exp(-z/a)
	else
		p <- (1-k*z/a)^(1/k)

	## k>0时，0 <= z <= a/k
	if(k > 0)
		p[z>a/k] <- 0

	return(p);
}


#' @title Calculate whether the GPD fit is good enough
#' @description Goodness of fit test for the generalized pareto distribution
#'   (gpd) P-value of the null hypothesis that the data comes from (or can be
#'   modeled with) the fitted gpd. Small p-values indicate a bad fit.
#' @parm p  cdf values for data fitted to the gpd (from gpcdf)
#' @parm k  estimated shape parameter of the gpd (from gpfit)
#' @return Pcm  P-value using Cramer-von Mises statistic
#' @return Pad  P-value using Anderson-Darling statistic (this gives more weight
#'   to observations in the tail of the distribution)
#' @return W2  Cramer-von Mises statistic
#' @return A2  Anderson-Darling statistic
#' @importFrom readr read_tsv
gpdgoft <- function(p, k){
	p <- sort(p)
	n <- length(p)
	i <- c(1:n)

	## Cramer-von Mises statistic
	W2 <- sum((p-((2*i-1)/(2*n)))^2) + 1/(12*n)
	## Anderson Darling statistic
	A2 <- -n - (1/n) * (sum((2*i-1) * matrix(log(p) + log(1-p[n+1-i]))))

	W2table <- c(0.046,0.067,0.094,0.115,0.136,0.165,0.187,0.239,
	             0.049,0.072,0.101,0.124,0.147,0.179,0.204,0.264,
	             0.053,0.078,0.111,0.137,0.164,0.2,0.228,0.294,
	             0.055,0.081,0.116,0.144,0.172,0.21,0.24,0.31,
	             0.057,0.086,0.124,0.153,0.183,0.224,0.255,0.33,
	             0.059,0.089,0.129,0.16,0.192,0.236,0.27,0.351,
	             0.062,0.094,0.137,0.171,0.206,0.254,0.291,0.38,
	             0.065,0.1,0.147,0.184,0.223,0.276,0.317,0.415,
	             0.069,0.107,0.159,0.201,0.244,0.303,0.349,0.458,
	             0.074,0.116,0.174,0.222,0.271,0.338,0.39,0.513)
	W2table <- matrix(W2table, ncol=8, nrow=10, byrow=TRUE)
	A2table <-c(0.339,0.471,0.641,0.771,0.905,1.086,1.226,1.559,
	            0.356,0.499,0.685,0.83,0.978,1.18,1.336,1.707,
	            0.376,0.534,0.741,0.903,1.069,1.296,1.471,1.893,
	            0.386,0.55,0.766,0.935,1.11,1.348,1.532,1.966,
	            0.397,0.569,0.796,0.974,1.158,1.409,1.603,2.064,
	            0.41,0.591,0.831,1.02,1.215,1.481,1.687,2.176,
	            0.426,0.617,0.873,1.074,1.283,1.567,1.788,2.314,
	            0.445,0.649,0.924,1.14,1.365,1.672,1.909,2.475,
	            0.468,0.688,0.985,1.221,1.465,1.799,2.058,2.674,
	            0.496,0.735,1.061,1.321,1.59,1.958,2.243,2.922)
	A2table <-matrix(A2table, ncol=8, nrow=10, byrow=TRUE)
	ktable <- matrix(c(0.9, 0.5, 0.2, 0.1, 0, -0.1, -0.2, -0.3, -0.4, -0.5), ncol = 1)
#	W2table <- data.frame(read_tsv(file = "W2table", col_names = FALSE))
#	A2table <- data.frame(read_tsv(file = "A2table", col_names = FALSE))
	ptable <- matrix(c(0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.001), ncol = 1)
	k <- max(0.5, k)

#	t <- approx(my.fun(ktable, W2table, k), ptable, W2, rule=2)$y
	Pcm <- max(min(approx(my.fun(ktable, W2table, k), ptable, W2, rule=2)$y, 1), 0)
#	t <- approx(my.fun(ktable, A2table, k), ptable, A2, rule=2)$y
	Pad <- max(min(approx(my.fun(ktable, A2table, k), ptable, A2, rule=2)$y, 1), 0)
#	Pcm <- max(min(approx(my.fun(ktable, W2table, k), ptable, W2)$y, 1), 0)
#	Pad <- max(min(approx(my.fun(ktable, A2table, k), ptable, A2)$y, 1), 0)

	return(list(Pcm = Pcm, Pad = Pad, W2 = W2, A2 = A2))
}

my.fun <- function(ktable, table, k){
	data <- c()
	for(i in 1:dim(table)[2])
		data <- c(data, approx(ktable, table[,i], xout=k, rule=2)$y)
	return(data)
}


