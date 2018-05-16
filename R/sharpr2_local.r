sharpr2_local <- function(data, s_a = 1000, s_m = 1, m_a = 0, len = 1, verbose = FALSE, mse = FALSE)
{
	nr <- nrow(data)

	max_b <- max(data$end) + 1

	T <- matrix(0, ncol = max_b, nrow = nr)
	for(i in 1:nr)
	{
		T[i,(data$start[i]:data$end[i])+1] <- 1
	}

	if(nr>1)
	{
		# L <- diag(data$length)
		L <- diag(len, nr)
	}else{
		L <- diag(1)
		L[1,1] <- len
	}

	# mu_a <- rep(m_a,max_b)
	# sigma_a <- diag(s_a,max_b)
	# sigma_a <- s_a*Diagonal(max_b)

	# sigma_m <- diag(rep(s_m,nr))
	
	L_i <- solve(L)
	
	if(ncol(T)>5000)
	{
		T <- Matrix(T, sparse=TRUE)
	}
			
	X <- L_i%*%T
	
	#cov_am <- sigma_a%*%t(T)%*%L_i
	#cov_m <- sigma_m + X%*%cov_am
	#cov_m_i <- solve(cov_m)
	
	lambda <- s_m/s_a

	#mu_m <- rep(m_a,nr)

	#est_a <- mu_a + cov_am%*%cov_m_i%*%(data$val-mu_m)
	
	s_d <- NA
	
	tryCatch({
	s_d <- svd(X)
	}, error = function(ex){}
	)
	
	if(is.na(s_d[1]))
	{
		if(verbose==TRUE)
		{warning("A tiled region is skipped due to an error in SVD.")}
		return(invisible(list(est_a = NA, var_nb = NA, mse = NA, lambda = NA, trh = NA)))
	}
	
	
	R <- s_d$u%*%diag(s_d$d)
	w1x <- s_d$v%*%solve(t(R)%*%R+lambda*diag(1,ncol(R)))%*%t(R)
	est_a <- w1x%*%data[,'val']
	
	
	sd_e <- NA
	#if(ci==TRUE)
	#{	
		
		# xtx <- t(X)%*%X
		# w_1 <- solve(xtx+lambda*diag(1,ncol(X)))
		# w1x <- w_1%*%t(X)
		h <- as.matrix(X%*%w1x)
		
		yxb <- as.matrix(data[,'val']-X%*%est_a)
		
		# sigma_e <- t(yxb)%*%yxb/(nrow(X)-2*sum(diag(h))+sum(diag(h%*%t(h))))
		sigma_e <- crossprod(yxb)/(nrow(X)-2*sum(diag(h))+sum(diag(tcrossprod(h))))
		# sd_nb <- sigma_e[1,1]*w1x%*%t(w1x)
		sd_nb <- sigma_e[1,1]*tcrossprod(w1x)
		
		if(mse==TRUE)
		{
			if(ncol(X)>5000)
			{
				X <- Matrix(X, sparse=TRUE)
			}
			## a little slow
			w <- w1x%*%X
			## a little slow
			w_b <- (w - diag(1,nrow(w)))%*%est_a
			# sd_e <- sqrt(diag(sigma_e[1,1]*w%*%t(w_1) + w_b%*%t(w_b)))
			## a little slow
			sd_e <- sqrt(diag(sd_nb + w_b%*%t(w_b)))
		}
	#}

	if(mse==TRUE)
	{
		return(invisible(list(est_a = est_a, mse = sd_e, var_nb = diag(sd_nb), lambda = lambda, trh = sum(diag(h)))))
	}else{
		return(invisible(list(est_a = est_a, mse = NA, var_nb = diag(sd_nb), lambda = lambda, trh = sum(diag(h)))))
	}
	
}
