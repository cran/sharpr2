call_sig_reg <- function(res, nr, threshold = 3.5, win = 10)
{
	if(class(res)!="sharpr2")
	{
		stop('The first argument must be an object obtained from sharpr2.')
	}
	
	re <- res$score[[nr]]
	
	sig_reg <- which(re$est_a - threshold*sqrt(re$var_nb) > 0)
	sig_reg_l <- length(sig_reg)
	if(sig_reg_l>1)
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
		sig_reg <- NA
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
		o_start <- as.numeric(strsplit(as.character(res$region[[nr]]),'-')[[1]][1])
		motif <- c(o_start + sig_reg[max_bp] - 1 -10, o_start + sig_reg[max_bp] -1 +10)
		sig_reg <- sig_reg + o_start - 1
	}
	
	res <- list(sig_reg=sig_reg, motif=motif)
	return(res)
}