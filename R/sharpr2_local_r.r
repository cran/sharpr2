sharpr2_local_r <- function(data, len = 1, verbose = FALSE, mse = FALSE, max_t = 1)
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
	
	L_i <- solve(L)
	
	if(ncol(T)>5000)
	{
		T <- Matrix(T, sparse=TRUE)
		X <- L_i%*%T
		X <- Matrix(X, sparse = TRUE)
	}else{
		X <- L_i%*%T
	}
	
	s_d <- NA
	
	tryCatch({
	s_d <- svd(X)
	}, error = function(ex){}
	)
	
	if(is.na(s_d[1]))
	{
		if(verbose==TRUE)
		{warning("A tiled region is skipped due to an error in SVD.")}
		return(invisible(list(est_a = NA, mse = NA, var_nb = NA, lambda = NA, trh = NA)))
	}
	
	min_r_v <- c()
	k_r_v <- c()
	# max_pc <- min(nrow(X),ncol(X))
	# max_pc <- sum(s_d$d>1e-10)
	
	max_pc <- max(1,ceiling(sum(s_d$d>1e-10)*max_t))
	e_v_2 <- s_d$d[1:max_pc]^2
	
	
	for(num_pc in 1:max_pc)
	{
	
		Q <- as.matrix(s_d$v[,1:num_pc])
		D <- s_d$d[1:num_pc]^2
		# alpha_r <- diag(1/D,num_pc)%*%t(Q)%*%t(X)%*%data[,'val']
		# res <- data[,'val'] - X%*%Q%*%alpha_r
		XQ <- as.matrix(X%*%Q)
		alpha_r <- diag(1/D,num_pc)%*%t(XQ)%*%data[,'val']
		res <- data[,'val'] - XQ%*%alpha_r
		sigma_a_e <-  (t(res)%*%res/(nrow(X)-num_pc))[1,1]
		# k_r <- num_pc*sigma_a_e/(t(alpha_r)%*%alpha_r)
		k_r <- num_pc*sigma_a_e/crossprod(alpha_r)
		min_r <- num_pc - sum((e_v_2^2)/((e_v_2+as.numeric(k_r))^2))
		min_r_v <- c(min_r_v,min_r)
		k_r_v <- c(k_r_v,k_r[1,1])
	}
	
	min_r_vs <- sort(min_r_v)
	ind_r <- 1
	num_pc <- match(min_r_vs[ind_r],min_r_v)
	lambda <- k_r_v[num_pc]
	
	## For some regions, lambda may be very small, which often happens when two problematic reads that are located nearly identically have opposite large values.
	## In this case, the last few eigenvalues may result in a very small sigma_a_e, and consequently a very small lambda that leads to unstable estimates. 
	## The following is used to give some correction for this.
	
	if((lambda<1e-05)&(verbose==TRUE))
	{warning("The data may contain two or more problematic reads located nearly identically that have opposite large values.")}
	
	while((max_pc<=max_b)&(lambda<1e-05))
	{
		ind_r <- ind_r + 1
		if(ind_r <= max_pc)
		{
			num_pc <- match(min_r_vs[ind_r],min_r_v)
			if(is.finite(k_r_v[num_pc]))
			{
				lambda <- k_r_v[num_pc]
			}else{
				lambda <- 0.0033
			}
		}else{
			lambda <- 0.0033
		}
		
	}
	 
	R <- s_d$u%*%diag(s_d$d)
	w1x <- s_d$v%*%solve(t(R)%*%R+lambda*diag(1,ncol(R)))%*%t(R)
	est_a <- w1x%*%data[,'val']
	
	sd_e <- NA
	sd_nb <- NA
	#if(ci==TRUE)
	#{
		# xtx <- t(X)%*%X
		# w_1 <- solve(xtx+lambda*diag(1,ncol(X)))
		# w1x <- w_1%*%t(X)
		h <- as.matrix(X%*%w1x)
		
		yxb <- as.matrix(data[,'val']-X%*%est_a)
		# sigma_e <- t(yxb)%*%yxb/(nrow(X)-sum(diag(h)))
		# sigma_e <- t(yxb)%*%yxb/(nrow(X)-1.25*sum(diag(h))+0.5)
		# sigma_e <- t(yxb)%*%yxb/(nrow(X)-2*sum(diag(h))+sum(diag(h%*%t(h))))
		sigma_e <- crossprod(yxb)/(nrow(X)-2*sum(diag(h))+sum(diag(tcrossprod(h))))
		# sd_nb <- sigma_e[1,1]*w1x%*%t(w1x)
		sd_nb <- sigma_e[1,1]*tcrossprod(w1x)
		
		if(mse==TRUE)
		{
			#if(ncol(X)>5000)
			#{
			#	X <- Matrix(X, sparse=TRUE)
			#}
	
			w <- as.matrix(w1x%*%X)
			w_b <- (w - diag(1,nrow(w)))%*%est_a
			
			# sd_e <- sqrt(diag(sigma_e[1,1]*w%*%t(w_1) + w_b%*%t(w_b)))
			sd_e <- sqrt(diag(sd_nb + w_b%*%t(w_b)))
		}
	#}
	if(mse==TRUE)
	{
		return(invisible(list(est_a = est_a, mse = sd_e, var_nb = sd_nb, lambda = lambda, trh = sum(diag(h)))))
	}else{
		return(invisible(list(est_a = est_a, mse = NA, var_nb = sd_nb, lambda = lambda, trh = sum(diag(h)))))
	}
	
}
