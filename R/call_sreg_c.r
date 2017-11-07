call_sreg_c <- function(re, alpha = 0.05, win = 10, o_start = 0)
{
	if(is.na(re$est_a[1]))
	{
		res <- list(sig_reg=NA, motif=NA, thres=NA)
		return(res)
	}
	
	ind <- abs(c(1,diff(re$var_nb[1,])))>1e-10
	est_sd <- sqrt(diag(re$var_nb)[ind])
	est <- re$est_a[ind]/est_sd
	l_reg <- length(est)
	if(l_reg<2)
	{
		res <- list(sig_reg=NA, motif=NA, thres=NA)
		return(res)
	}
	
	scale_m <- est_sd%*%t(est_sd)
	corr_m <- re$var_nb[ind,ind]/scale_m
	
	thres <- NA
	ntry <- 0
	while((is.na(thres))&(ntry<5))
	{
		tryCatch({
			if(l_reg<=1000)
			{
				thres <- qmvnorm(1-alpha,interval=c(-10, 10+2*ntry),tail = c("lower.tail"), mean=rep(0,l_reg),corr=corr_m[1:l_reg,1:l_reg])$quantile
			}else{
				thres <- qmvnorm(1-alpha,interval=c(-10, 10+2*ntry),tail = c("lower.tail"), mean=rep(0,1000),corr=corr_m[1:1000,1:1000])$quantile
			}
		}, error = function(ex){}
		)
		ntry <- ntry + 1
	}
	
	if(is.na(thres)|(thres<0))
	{
		res <- list(sig_reg=NA, motif=NA, thres=NA)
		return(res)
	}
	
	sig_reg <- which(re$est_a/sqrt(diag(re$var_nb)) - thres > 0)
	sig_reg_l <- length(sig_reg)
	if((sig_reg_l>1)&(win>1))
	{
		keep <- rep(TRUE, sig_reg_l)
		p <- 1
		count <- 1
		for(cum in 2:sig_reg_l)
		{
			if(sig_reg[cum]==(sig_reg[cum-1]+1))
			{
				count <- count + 1
				if(cum==sig_reg_l)
				{
					if(count<=win)
					{
						keep[p:cum] <- FALSE
					}
				}
			}else{
				if(count<=win)
				{
					keep[p:(cum-1)] <- FALSE
				}
				count <- 1
				p <- cum 
			}
			
		}
		
		sig_reg <- sig_reg[keep]
		if(length(sig_reg)<=win)
		{
			sig_reg <- NA
		}
	}else{
		if(sig_reg_l<1)
		{sig_reg <- NA}
	}
	
	motif <- NA
	
	if(!is.na(sig_reg[1]))
	{
		max_mean <- c()
		for(bp in 1:length(sig_reg))
		{
			max_mean <- c(max_mean,mean(re$est_a[max(1,sig_reg[bp]-10):min(sig_reg[bp]+10,length(re$est_a)),1]))
		}
		max_bp <- which(max_mean==max(max_mean))
		max_bp <- max_bp[ceiling(length(max_bp)/2)]
		# o_start <- as.numeric(strsplit(as.character(res$region[[nr]]),'-')[[1]][1])
		motif <- c(o_start + sig_reg[max_bp] - 1 -10, o_start + sig_reg[max_bp] - 1 +10)
		sig_reg <- sig_reg + o_start - 1
	}
	
	res <- list(sig_reg=sig_reg, motif=motif, thres = thres)
	return(res)
}