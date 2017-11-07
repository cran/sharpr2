#' plot.sharpr2
#'
#' For an ATAC-STAR dataset on a given chromosome, this function plots the estimated scores (with s.e. if available) for a tile region.
#' @param re A object returned from the sharpr2 funciton. 
#' @param tr A integer indicating which region in re to be plotted.
#' @param loess An indicator for whether the loess method is used for smoothing in plotting the scores from sharpr2. The standard errors are not plotted when loess is used.
#' @param add An indicator for whether to add the new plot to the existing one.
#' @param xlab The label for the x-axis. The default is 'Position'.
#' @param ylab The label for the y-axis. The default is 'Regulatory Score'.
#' @param cicol The color for CIs. The default is 'orange'.
#' @param cimcol The color for filling the regions within CIs. The default is 'grey'.
#' @param sreg An indicator for whether to highlight the significant regions. The default is TRUE.
#' @param ... Other parameters for plot.
#' @keywords sharpr2 ATAC-STAR
#' @export
#' @examples
#' # plot(re,1,loess=FALSE)


plot.sharpr2 <- function(x, tr, unc = 'CI', loess=FALSE, add=FALSE, xlab='Position', ylab='Regulatory Score', cicol = "orange", cimcol = 'grey', sreg = TRUE, ...)
{
	if(is.na(tr))
	{
		stop("The argument tr must be an integer.")
	}
	
	region <- as.numeric(unlist(strsplit(x$region[[tr]],'-')))
	
	if(is.na(x$score[[tr]]$est_a[1]))
	{
		stop("There is no result for this tiled region.")
	}
	
	if(var(x$score[[tr]]$est_a)==0)
	{
		if(add==FALSE)
		{
			plot.default(region[1]:region[2],x$score[[tr]]$est_a,type='l',xlab=xlab, ylab=ylab, ...)
		}else{
			lines(region[1]:region[2],x$score[[tr]]$est_a,xlab=xlab, ylab=ylab, ...)
		}
		
	}else{
		if(loess==TRUE)
		{	
			if(add==FALSE)
			{
				plot.default(loess.smooth(region[1]:region[2],x$score[[tr]]$est_a,degree = 2),type='l',xlab=xlab, ylab=ylab, ...)
			}else{
				lines(loess.smooth(region[1]:region[2],x$score[[tr]]$est_a,degree = 2),xlab=xlab, ylab=ylab, ...)
			}
		}else{
			if(add==FALSE)
			{
				if(is.na(x$score[[tr]]$var_nb[1]))
				{
					plot.default(region[1]:region[2],x$score[[tr]]$est_a,type='l',xlab=xlab, ylab=ylab, ...)
				}else{
					if(unc=='MSE')
					{
						if(!is.na(x$score[[tr]]$mse[1]))
						{
							lower <- x$score[[tr]]$est_a - x$score[[tr]]$mse
							upper <- x$score[[tr]]$est_a + x$score[[tr]]$mse
						}else{
							warning("MSEs are not given in the results.")
						}
						
					}else{
						if(unc=='CI')
						{
							lower <- x$score[[tr]]$est_a - 1.96*sqrt(x$score[[tr]]$var_nb)
							upper <- x$score[[tr]]$est_a + 1.96*sqrt(x$score[[tr]]$var_nb)
						}	
					}		
					
					plot.default(region[1]:region[2],x$score[[tr]]$est_a,type='l', ylim=c(min(lower),max(upper)), xlab=xlab, ylab=ylab, ...)
					
					lines(region[1]:region[2], lower, col = cicol ,lty = 2 , lwd = 0.6)
					lines(region[1]:region[2], upper, col = cicol ,lty = 2 , lwd = 0.6)
					polygon(c(region[1]:region[2], rev(region[1]:region[2])),c(upper, rev(lower)),col=cimcol,border = NA, lty=3, density=20)
					lines(region[1]:region[2],x$cutoff[[tr]]*sqrt(x$score[[tr]]$var_nb), col=rgb(0,1,0,alpha=0.5))
					if(!is.na(x$sig_reg[[tr]][1]))
					{
						if(sreg==TRUE)
						{
							points(x=x$sig_reg[[tr]], y=rep(0,length(x$sig_reg[[tr]])),col='red',cex=1,pch=20)
							n_sreg <- split(x$sig_reg[[tr]], cumsum(c(1, diff(x$sig_reg[[tr]]) != 1)))
							for(j in 1:length(n_sreg))
							{polygon(c(n_sreg[[j]][1],n_sreg[[j]][length(n_sreg[[j]])],n_sreg[[j]][length(n_sreg[[j]])],n_sreg[[j]][1]), c(min(lower),min(lower),max(upper),max(upper)), col=rgb(1, 0, 0,0.4), border='black', lty = c("dashed"))}
							#points(x=re$motif[[tr]][1]:re$motif[[tr]][2], y=rep(0,re$motif[[tr]][2]-re$motif[[tr]][1]+1),col='blue',cex=1,pch=20)
							polygon(c(x$motif[[tr]],rev(x$motif[[tr]])), c(min(lower),min(lower),max(upper),max(upper)), col=rgb(0, 0, 1,0.2), border=NA)
							# lines(region[1]:region[2],re$cutoff[[tr]]*sqrt(re$score[[tr]]$var_nb), col=rgb(0,1,0,alpha=0.5))
						}
					}
				}
			}else{
				if(is.na(x$score[[tr]]$var_nb[1]))
				{
					lines(region[1]:region[2],x$score[[tr]]$est_a,xlab=xlab, ylab=ylab, ...)
				}else{
					if(unc=='MSE')
					{
						if(!is.na(x$score[[tr]]$mse[1]))
						{
							lower <- x$score[[tr]]$est_a - x$score[[tr]]$mse
							upper <- x$score[[tr]]$est_a + x$score[[tr]]$mse
						}else{
							warning("MSEs are not given in the results.")
						}
					}else{
						if(unc=='CI')
						{
							lower <- x$score[[tr]]$est_a - 1.96*sqrt(x$score[[tr]]$var_nb)
							upper <- x$score[[tr]]$est_a + 1.96*sqrt(x$score[[tr]]$var_nb)
						}
					}
					
					lines(region[1]:region[2],x$score[[tr]]$est_a,type='l', ylim=c(min(lower),max(upper)), xlab=xlab, ylab=ylab, ...)
					
					lines(region[1]:region[2], lower, col = cicol ,lty = 2 , lwd = 0.6)
					lines(region[1]:region[2], upper, col = cicol ,lty = 2 , lwd = 0.6)
					polygon(c(region[1]:region[2], rev(region[1]:region[2])),c(upper, rev(lower)),col=cimcol,border = NA, lty=3, density=20)
					lines(region[1]:region[2], x$cutoff[[tr]]*sqrt(x$score[[tr]]$var_nb), col=rgb(0,1,0,alpha=0.5))
					
					if(!is.na(x$sig_reg[[tr]][1]))
					{
						if(sreg==TRUE)
						{
							points(x=x$sig_reg[[tr]], y=rep(0,length(x$sig_reg[[tr]])),col='red',cex=1,pch=20)
							n_sreg <- split(x$sig_reg[[tr]], cumsum(c(1, diff(x$sig_reg[[tr]]) != 1)))
							for(j in 1:length(n_sreg))
							{polygon(c(n_sreg[[j]][1],n_sreg[[j]][length(n_sreg[[j]])],n_sreg[[j]][length(n_sreg[[j]])],n_sreg[[j]][1]), c(min(lower),min(lower),max(upper),max(upper)), col=rgb(1, 0, 0,0.4), border='black', lty = c("dashed"))}
							#points(x=re$motif[[tr]][1]:re$motif[[tr]][2], y=rep(0,re$motif[[tr]][2]-re$motif[[tr]][1]+1),col='blue',cex=1,pch=20)
							polygon(c(x$motif[[tr]],rev(x$motif[[tr]])), c(min(lower),min(lower),max(upper),max(upper)), col=rgb(0, 0, 1,0.2), border=NA)
							# lines(region[1]:region[2], re$cutoff[[tr]]*sqrt(re$score[[tr]]$var_nb), col=rgb(0,1,0,alpha=0.5))
						}
					}
				
				}
			}
			
		}
	}
}
