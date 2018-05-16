#' sharpr2
#'
#' For a HiDRA dataset on a given chromosome, this function calls tiled regions (the regions covered by at least one read), and calculates regulatory scores for each tiled region. The regulatory scores are based on standardized log(RNA/PLASMID).
#' @param data A data.frame for an HiDRA dataset for one chromosome. The data.frame must contain four columns: 'start', 'end', 'PLASMID', 'RNA' 
#' @param l_min the minimum length for filtering a read. The default is 150.  
#' @param l_max the maximum length for filtering a read. The default is 600.
#' @param f_rna the minimum RNA count for filtering a read. The default is 10.
#' @param f_dna the minimum PLASMID count for filtering a read. The default is 0.
#' @param s_a A variance hyperparameter in the prior for the latent regulatory scores. The default is 300.
#' @param verbose An indicator of whether to show processing information. The default is FALSE.
#' @param auto An indicator of whether to automatically estimate the hyperparameters from the data using a ridge regression. The default is TRUE.
#' @param sig An indicator of whether to identify significant regions for the estimated scores. Only valid if auto=TRUE. The default is TRUE.
#' @param len An indicator of whether to model log(RNA/PLASMID) of each read as the average or the sum of the latent regulatory scores. The default is FALSE, which is the sum.
#' @param alpha A threshold to call the significant region. The default is 0.05.
#' @param win A window size for removing sporadic significant regions. If a significant consecutive region is small than win, it will be treated as false signals. The default is 5.
#' @keywords sharpr2 HiDRA
#' @return score: the regulatory scores for each tiled region.
#' @return region: the start and end positions for each tiled region.
#' @return n_reg: total number of tiled regions.
#' @return n_read: the number of reads in each tiled region.
#' @export
#' @examples
#' # sharpr2(data)


sharpr2 <- function(data, l_min = 150, l_max = 600, f_rna = 10, f_dna = 0, s_a = 300, verbose = FALSE, auto = TRUE, sig = TRUE, len = FALSE, alpha = 0.05, win = 5, mse = FALSE, max_t = 1)
{
	data$length <- data$end - data$start + 1
	
	data_c <- data[which((data$length>l_min)&(data$length<l_max)),]
	data_c <- data_c[which((data_c$RNA>f_rna)),]
	data_c <- data_c[which((data_c$PLASMID>f_dna)),]
	
	if((f_dna<0)|(f_rna<0))
	{
		stop("The arguments f_dna and f_rna must be non-negative.")
	}
	
	if((alpha<=0)|(alpha>=1))
	{
		stop("The argument alpha must be between 0 and 1.")
	}
	
	if(max_t>1)
	{max_t <- 1}
	
	if(verbose == TRUE)
	{
		cat('Total reads after filtering: ', nrow(data_c), '\n')
	}
	
	# data_chr1 <- data_c[which(data_c$chr=='chr1'),]
	data_c <- data_c[order(data_c$start),]
	
	t_r <- call_tile_reg(data_c)
	
	n_s <- t_r$num_r
	reg_set <- t_r$tile_reg
	size <- t_r$size
	
	if(verbose == TRUE)
	{
		cat('Call tiled regions: ', n_s, '\n')
	}
	
	
	re_list <- vector('list',n_s)
	re_reg <- vector('list',n_s)
	sig_regs <- vector('list',n_s)
	motifs <- vector('list',n_s)
	cutoffs <- vector('list',n_s)
	data_c$val <- log(data_c$RNA/data_c$PLASMID)
	pnull <- unlist(t_r$tile_reg[t_r$size<5])
	sd_pn <- sd(data_c$val[pnull])
	if((length(pnull)>2)&(!is.na(sd_pn))&(sd_pn>0))
	{
		t_mean <- mean(data_c$val[pnull])
		t_sd <- sd(data_c$val[pnull])
	}else{
		t_mean <- mean(data_c$val)
		t_sd <- sd(data_c$val)
	}
	if((t_sd==0)&(is.nan(t_sd)))
	{	
		warning("There is no variation of the values in the library or the data contains Inf.")
		t_sd <- 1
	}
	dna_mean <- mean(data_c$PLASMID)
	len_mean <- mean(data_c$length)
	
	if(verbose == TRUE)
	{
		cat('Calculating regulatory scores ... ', '\n')
	}
	
	if(n_s==0)
	{
		stop("No reads pass the quality control.")
	}
	
	for(i in 1:n_s)
	{
		id <- reg_set[i]
		t_d <- data_c[id[[1]],]
		re_reg[[i]] <- paste(min(t_d$start),'-',max(t_d$end),sep='')
		
		o_start <- t_d$start[1]
		## shift to start=0
		t_d$end <- t_d$end - t_d$start[1]
		t_d$start <- t_d$start - t_d$start[1]
		#t_d$PLASMID <- apply(t_d[,2:6],1,mean)
		t_d$PLASMID <- ifelse(t_d$PLASMID==0,1,t_d$PLASMID)
		#t_d$RNA <- apply(t_d[,7:11],1,mean)
		# t_d$val <- log(t_d$RNA/t_d$PLASMID)
		t_d$val <- (log(t_d$RNA/t_d$PLASMID)-t_mean)/t_sd
		t_d$w <- t_d$PLASMID/dna_mean
		
		if(size[i]>1)
		{
			if(len == FALSE)
			{
				len_t <- len_mean
			}else{
				len_t <- t_d$length
			}
			
			if(auto==FALSE)
			{
				re <- sharpr2_local(t_d[,c('start','end','length','val','w')], s_m=1, s_a=s_a, m_a=0, len = len_t, mse = mse)
				sig_reg <- NA
				motif <- NA
				thres <- NA
			}else{
				re <- sharpr2_local_r(t_d[,c('start','end','length','val','w')], len = len_t, verbose = verbose, mse = mse, max_t = max_t)
				# re <- sharpr2_local(t_d[,c('start','end','length','val','w')], s_m=1, s_a=1000, m_a=0, weight=FALSE, ci=ci)
				if(!is.na(re$est_a[1]))
				{
					if(sig==TRUE)
					{
						res <- call_sreg_c(re, alpha, win, o_start)
						sig_reg <- res$sig_reg
						motif <- res$motif
						thres <- res$thres
					}else{
						sig_reg <- NA
						motif <- NA
						thres <- NA
					}
					re$var_nb <- diag(re$var_nb)
				}else{
					sig_reg <- NA
					motif <- NA
					thres <- NA
				}
			}
			
			if(verbose == TRUE)
			{
				print(i)
			}
		}else{
			re <- list(est_a = rep(t_d$val,t_d$end-t_d$start+1), mse = NA, var_nb = NA, lambda = NA)
			sig_reg <- NA
			motif <- NA
			thres <- NA
		}
		
		re_list[[i]] <- re
		sig_regs[[i]] <- sig_reg
		motifs[[i]] <- motif
		cutoffs[[i]] <- thres
	}
	
	if(verbose == TRUE)
	{
		cat('Done!', '\n')
	}
	
	re <- list(score = re_list, region = re_reg, n_reg = n_s, n_read = size, sig_reg = sig_regs, motif = motifs, cutoff = cutoffs)
	class(re) <- "sharpr2"
	return(invisible(re))
}

