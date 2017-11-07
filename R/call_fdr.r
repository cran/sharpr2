#' call_fdr
#'
#' 
#' @param whole_re a list of objects of class sharpr2 for each chromosome.  
#' @param thres_tr the threshold for the size of tiled reigons used for calculate FDR-adjusted p-values. The default value is 10.  
#' @param method the method for calculating FDR-adjusted p-values. See the function 'p.adjust' for more details abount the method. The default is 'BH'.
#' @keywords sharpr2 HiDRA
#' @return gfdr: a result table of FDR-adjusted p-values.
#' @export
#' @examples
#' # call_fdr(data)


call_fdr <- function(whole_re, thres_tr = 10, method = 'BH')
{
	num_chr <- length(whole_re)
	
	if(num_chr==0)
	{
		stop('No chromosome is found in the data.')
	}
	
	p_v <- c()
	i_reg_s <- c()
	i_reg_e <- c()
	num_tr <- c()
	chr <- c()
	n_read <- c()
	for(i in 1:num_chr)
	{
		# chr_r <- alpha005_win1_res[[i]]
		chr_r <- whole_re[[i]]
		all_r <- strsplit(as.character(chr_r$region),'-')
		lr <- which(chr_r$n_read>thres_tr)
		for(j in lr)
		{
			re <- chr_r$score[[j]]
			reg_t <- as.numeric(all_r[[j]])
			o_start <- reg_t[1]
			ind <- abs(c(1,diff(re$var_nb)))>1e-10
			# est <- pnorm(re$est_a[ind]/sqrt(re$var_nb[ind]),lower.tail=FALSE)
			est <- pt(re$est_a[ind]/sqrt(re$var_nb[ind]),df=chr_r$n_read[j]-re$trh,lower.tail=FALSE)
			p_v <- c(p_v, est)
			s_t <- o_start + which(ind) - 1
			i_reg_s <- c(i_reg_s, s_t)
			if(length(s_t)<2)
			{
				i_reg_e <- c(i_reg_e, reg_t[2]) 
			}else{
				i_reg_e <- c(i_reg_e, c(s_t[2:length(s_t)]-1, reg_t[2]))
			}
			num_tr <- c(num_tr, rep(j, length(est)))
			chr <- c(chr, rep(i, length(est)))
			n_read <- c(n_read, rep(chr_r$n_read[[j]], length(est)))
		}
	}

	gfdr <- data.frame(p=p_v, start=i_reg_s, end=i_reg_e, ntr=num_tr, chr=chr, size_tr=n_read)
	gfdr$fdr <- p.adjust(gfdr$p, method = method)
	
	return(gfdr)

}
