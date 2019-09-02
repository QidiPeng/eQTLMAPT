REALMIN <- 2.225073858507201e-308
EPS <- 2.220446049250313e-16

## Implement matlab eps function
eps <- function(x){
	if(abs(x) < REALMIN)
		return(2^(-1074))

	E <- log(abs(x), base=2)
	return(2^(E-52))
}


#' @title Fitting the generalized pareto distribution using 'ML'
#' @parm z  exceedances
#' @return Pharmat  estimated shape and scale parameter
gpdfit <- function(z){
	TolBnd <- 1e-10
	TolFun <- 1e-10
	TolX <- 1e-10

	parmhat = c(1, 0)
	res <- SimplexMethod(parmhat, MaxFunEvals = 2500, MaxIter = 1000, TolFun = TolFun, TolX = TolX, z)
	parmhat <- res$x
	err <- res$exitflag

	c <- 0
	while(err == 0 & c < 10){
		TolBnd <- TolBnd * 10
		TolFun <- TolFun * 10
		TolX <- TolX * 10
		res <- SimplexMethod(parmhat, MaxFunEvals = 2500, MaxIter = 1000, TolFun = TolFun, TolX = TolX, z)
		parmhat <- res$x
		err <- res$exitflag
		c <- c + 1
	}

	## k==1 (boundary maximum)
	if(abs(parmhat[2] - 1) < TolBnd){
		parmhat <- c(NaN, NaN)
	}

	return(parmhat)
}



#' @title Multidimensional unconstrained nonlinear minimization (Nelder-Mead)
#' @parm x  Initial point.
#' @parm MaxFunEvals  Maximum objective function calculations.
#' @parm MaxIter  The maximum number of iterations.
#' @parm TolFun  The termination tolerance for function values.
#' @parm TolX  The termination tolerance for the current point x.
#' @parm varargin  Additional parameters of the nll function (non-optimized
#'   parameters).
#' @return x  Minimum objective function value point.
#' @return fval  Minimum objective function value.
#' @return exitflag  Return flag. If the result is found, it returns 1; if the
#'   number of function calculations or the number of iterations reaches the set
#'   value, there is still no result, and 0 is returned.
SimplexMethod <- function(x, MaxFunEvals = 'default', MaxIter = 'default', TolFun = 1e-10, TolX = 1e-10, varargin){
	options(digits = 22)
	n <- length(x)
	if(MaxFunEvals == 'default')
		MaxFunEvals = 200 * n
	if(MaxIter == 'default')
		MaxIter = 200 * n

	## Initialize parameters
	rho <- 1
	chi <- 2
	psi <- 0.5
	sigma <- 0.5
	onesn <- matrix(data=1, nrow=1, ncol=n)
	two2np1 <- c(2:n+1)
	one2n <- c(1:n)

	## Set up a simplex near the initial guess
	xin <- matrix(x)
	v <- matrix(data=0, nrow=n, ncol=n+1)
	fv <- matrix(data=0, nrow=1, ncol=n+1)
	v[,1] <- xin
##	x <- xin[,1]
	fv[1,1] <- nll(x, varargin)
	func_evals <- 1
	itercount <- 0

	usual_delta <- 0.05
	zero_term_delta <- 0.00025
	for(j in 1:n){
		y <- xin
		if(y[j] == 0)
			y[j] = zero_term_delta
		else
			y[j] = (1 + usual_delta) * y[j]

		v[,j+1] = y
		x <- y[,1]
		f <- nll(x, varargin)
		fv[1,j+1] <- f
	}
	idx <- order(fv)
	fv <- fv[,idx]
	v <- v[,idx]

	itercount <- itercount + 1
	func_evals <- n + 1

	exitflag <- 1
	## iteration
	while(func_evals < MaxFunEvals && itercount < MaxIter){
		if((max(abs(fv[1]-fv[two2np1])) <= max(TolFun,10*eps(fv[1]))) &&
				(max(abs(v[,two2np1]-v[,onesn])) <= max(TolX,10*eps(v[,1]))))
			break;

		##Compute the reflection point
		xbar <- rowSums(v[,one2n]) / n
		xr <- 2 * xbar - v[,n+1]
		fxr <- nll(xr,varargin)
		func_evals <- func_evals + 1

		if(fxr < fv[1]){
			##Calculate the expansion point
			xs <- (1 + chi) * xbar - chi * v[,n+1]
			fxs <- nll(xs, varargin)
			func_evals <- func_evals + 1
			if(fxs < fxr){
				v[,n+1] <- xs
				fv[n+1] <- fxs
				how <- 'expand'
			}else{
				v[,n+1] <- xr
				fv[n+1] <- fxr
				how <- "refelect"
			}
		}else{
			if(fxr < fv[n]){
				v[,n+1] <- xr
				fv[n+1] <- fxr
				how <- "refelect"
			}else{
				##Perform contraction
				if(fxr < fv[n+1]){
					##Perform an outside contraction
					xc <- (1 + psi) * xbar - psi * v[,n+1]
					fxc <- nll(xc, varargin)
					func_evals <- func_evals + 1
					if(fxc <= fxr){
						v[,n+1] <- xc
						fv[n+1] <- fxc
						how <- 'contract outside'
					}else{
						##perform a shrink
						how <- 'shrink'
					}
				}else{
					##Perform an inside contraction
					xcc <- psi * xbar + psi * v[,n+1]
					fxcc = nll(xcc, varargin)
					func_evals <- func_evals + 1
					if(fxcc < fv[n+1]){
						v[,n+1] <- xcc
						fv[n+1] <- fxcc
						how <- 'contract inside'
					}else{
						##perform a shrink
						how <- 'shrink'
					}
				}

				if(how == 'shrink'){
					for(j in two2np1){
						v[,j] <- v[,1] + sigma * (v[,j] - v[,1])
						fv[j] <- nll(v[,j], varargin)
					}
					func_evals <- func_evals + n
				}
			}
		}

		idx <- order(fv)
		fv <- fv[idx]
		v <- v[,idx]
		itercount <- itercount + 1
	}#end of while

	x <- v[,1]
	fval <- fv[1]
	if(func_evals >= MaxFunEvals || itercount >= MaxIter)
		exitflag <- 0

	return(list(x = x, fval = fval, exitflag = exitflag))
}


nll <- function(parmhat, z){
	a <- parmhat[1]
	k <- parmhat[2]
	n <- length(z)
	m <- max(z)

	if(k > 1)
		L <- -Inf
	else if(a < max(0,k*m))
		L <- -Inf
	else{
		if(abs(k) < EPS){ #k==0(<EPS)
			L = -n*log(a)-1/a*sum(z)
		}
		else{
			L = -n*log(a) + (1/k-1)*sum(log(1-k*z/a))
		}
	}
	L <- -L
	return(L)
}

